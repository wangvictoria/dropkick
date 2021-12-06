from django.shortcuts import render, get_object_or_404, redirect
from dropkickApp.models import MyFile, CustomParam
from django.views import generic
from .forms import UploadFileForm, CheckboxForm, CustomForm, ScoreForm
from django.http import HttpResponse, StreamingHttpResponse, FileResponse
from django.core.files.storage import FileSystemStorage
import csv
import os
import zipfile
from io import BytesIO
from django.core.exceptions import ValidationError
from django.contrib import messages

import scanpy as sc; sc.set_figure_params(color_map="viridis", frameon=False)
import dropkick as dk
import matplotlib.pyplot as plt; plt.switch_backend("Agg")
import io, base64, urllib
import numpy as np
import pandas as pd


def qc_plot(adata):
    # plot QC metrics
    adata = dk.recipe_dropkick(adata, n_hvgs=None, X_final="raw_counts")
    qc_plt = dk.qc_summary(adata)
    
    # display chart
    buf = io.BytesIO()
    qc_plt.savefig(buf, format = 'png')
    qc_plt.savefig('media/qc_plot.png')
    buf.seek(0)
    string = base64.b64encode(buf.read())
    uri = urllib.parse.quote(string)
    return uri

def labels(adata, min_genes, mito_names, n_ambient, n_hvgs, thresh_methods, alphas, max_iter, seed):
    adata_model = dk.dropkick(
        adata, 
        min_genes=min_genes, 
        mito_names=mito_names, 
        n_ambient=n_ambient,
        n_hvgs=n_hvgs,
        thresh_methods=thresh_methods,
        alphas=alphas,
        max_iter=max_iter,
        n_jobs=5,
        seed=seed)
    
    # display coefficient plot
    coef_plt = dk.coef_plot(adata)
    buf_coef = io.BytesIO()
    coef_plt.savefig(buf_coef, format = 'png')
    coef_plt.savefig('media/coef_plot.png')
    buf_coef.seek(0)
    string_coef = base64.b64encode(buf_coef.read())
    uri_coef = urllib.parse.quote(string_coef)
    
    # display score plot
    adata_score = dk.recipe_dropkick(adata, n_hvgs=None, verbose=False, filter=True, min_genes=50)
    score_plt = dk.score_plot(adata_score)
    buf_score = io.BytesIO()
    score_plt.savefig(buf_score, format = 'png')
    score_plt.savefig('media/score_plot.png')
    buf_score.seek(0)
    string_score = base64.b64encode(buf_score.read())
    uri_score = urllib.parse.quote(string_score)
    
    return adata, uri_score, uri_coef

def param_assignment(instance):
    if not instance.min_genes:
        instance.min_genes = 50
        instance.save()
    if not instance.mito_names:
        instance.mito_names = '^mt-|^MT-'
        instance.save()
    if not instance.n_ambient:
        instance.n_ambient = 10
        instance.save()
    if not instance.n_hvgs:
        instance.n_hvgs = 2000
        instance.save()
    if not instance.score_thresh:
        instance.score_thresh = 0.5
        instance.save()
    if not instance.alphas:
        instance.alphas = '0.1'
        instance.save()
    if not instance.max_iter:
        instance.max_iter = 2000
        instance.save()
    if not instance.seed:
        instance.seed = 18
        instance.save()

def index(request):
    """View function for home page of site."""
    form = CustomForm(request.POST or None, initial = {'min_genes': 50, 'mito_names': 'mt'})
    model = CustomParam
    # upload file
    if request.method == 'POST':
        if 'document' in request.FILES:
            if form.is_valid():
                form.save()
                instance = model.objects.last()
            
                uploaded_file = request.FILES['document']
                if uploaded_file.name.endswith('.csv'):
                    fs = FileSystemStorage()
                    fs.save(uploaded_file.name, uploaded_file)
                    os.rename('media/' + uploaded_file.name, 'media/sample.csv')
                    instance.csv_bool = True
                    instance.save()
                    param_assignment(instance)
                    return redirect(process)
                elif uploaded_file.name.endswith('.h5ad'):
                    fs = FileSystemStorage()
                    fs.save(uploaded_file.name, uploaded_file)
                    os.rename('media/' + uploaded_file.name, 'media/sample.h5ad')
                    instance.h5ad_bool = True
                    instance.save()
                    param_assignment(instance)
                    return redirect(process)
                elif uploaded_file.name.endswith('.tsv'):
                    fs = FileSystemStorage()
                    fs.save(uploaded_file.name, uploaded_file)
                    os.rename('media/' + uploaded_file.name, 'media/sample.tsv')
                    instance.tsv_bool = True
                    instance.save()
                    param_assignment(instance)
                    return redirect(process)
                else:
                    messages.error(request,'Please upload a file of CSV, H5AD, or TSV type')
                
            else:
                form = CustomForm(request.POST or None)
        else:
            messages.error(request,'Please select a file.')
#     else:
#         form = CheckboxForm()
        
        
    return render(request,'index.html', context = {
        'form': form
    })

def process(request):
    context = {
            'title': None, 'counts_text': None, 'counts_false': None, 'counts_true': None, 
            'qc_text': None, 'score_text': None, 'coef_text': None, 'labels_text': None,
            'qc_plot': None, 'score_plot': None, 'coef_plot': None, 'labels': None,
        }
    
    model = CustomParam
    instance = model.objects.last()
    if instance.csv_bool:
        adata = sc.read('media/sample.csv')
    if instance.h5ad_bool:
        adata = sc.read('media/sample.h5ad')
    if instance.tsv_bool:
        adata = sc.read('media/sample.tsv')
    
    # label data results
    context['title'] = 'Your Results'
    if request.method == 'POST':
        form = ScoreForm(request.POST or None)
        model = CustomParam
        if form.is_valid():
            instance.score_thresh = form.cleaned_data.get('score_thresh')
            instance.save()
            return redirect(calc_score_thresh)
                
        else:
            form = ScoreForm()
    if instance.qc_plot:
        # qc_plot checkbox was checked
        context['qc_text'] = 'QC Plot'
        context['qc_plot'] = qc_plot(adata)
        
    if instance.dropkick:
        # filter checkbox was checked

        # run dropkick
        context['counts_text'] = 'Droplets Inventory'
        context['score_text'] = 'Score Plot'
        context['coef_text'] = 'Coefficient Plot'
        context['labels_text'] = 'Dropkick Labels'
        
        alphas_list = instance.alphas.split(',')
        alphas = [float(x) for x in alphas_list]


        df, context['score_plot'], context['coef_plot'] = labels(
            adata, instance.min_genes, instance.mito_names, instance.n_ambient, instance.n_hvgs, instance.thresh_methods, alphas,
            instance.max_iter, instance.seed)
        
        context['score_thresh'] = instance.score_thresh

        df.obs['dropkick_label'] = df.obs['dropkick_score'] > instance.score_thresh

        context['counts_false'] = df.obs['dropkick_label'].value_counts()[0]
        context['counts_true'] = df.obs['dropkick_label'].value_counts()[1]

        # convert dataframe to csv
        df.obs.to_csv('media/dropkick_labels.csv')

        # convert to h5ad file
        adata.write('media/dropkick_filter.h5ad', compression='gzip')

        # output counts and genes matrices
        data_out = adata[df.obs['dropkick_label']==True]
        data = pd.DataFrame(data_out.X.toarray())
        data.to_csv('media/dropkick_counts.csv', header=False, index=False)
        pd.DataFrame(data_out.var_names).to_csv('media/dropkick_genes.csv', header=False, index=False)

    return render(request, 'process.html', context)


def calc_score_thresh(request):
    context = {
        'score_thresh': None, 'title': None, 'qc_text': None, 'counts_text': None, 'counts_false': None, 'counts_true': None,
    }
    form = ScoreForm(request.POST or None)
    model = CustomParam
    instance = model.objects.last()
    context['title'] = 'Your Results'
    if request.method == 'POST':
        if form.is_valid():
            instance.score_thresh = form.cleaned_data.get('score_thresh')
            instance.save()
            return redirect(calc_score_thresh)
                
        else:
            form = ScoreForm()
    if instance.qc_plot:
        # qc_plot checkbox was checked
        context['qc_text'] = 'QC Plot'
        
    if instance.dropkick:
        # filter checkbox was checked

        # run dropkick
        context['counts_text'] = 'Droplets Inventory'
        context['score_text'] = 'Score Plot'
        context['coef_text'] = 'Coefficient Plot'
        context['labels_text'] = 'Dropkick Labels'
    
        score_thresh = instance.score_thresh
        context['score_thresh'] = score_thresh
                
        df = sc.read('media/dropkick_filter.h5ad')
                
        df.obs['dropkick_label'] = df.obs['dropkick_score'] > score_thresh

        context['counts_false'] = df.obs['dropkick_label'].value_counts()[0]
        context['counts_true'] = df.obs['dropkick_label'].value_counts()[1]
        

        # convert dataframe to csv
        df.obs.to_csv('media/dropkick_labels.csv')

        # convert to h5ad file
        df.write('media/dropkick_filter.h5ad', compression='gzip')
    return render(request, 'score_thresh.html', context)

def download_csv(request):
    file = open('media/dropkick_labels.csv', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)

    # decide the file name
    response['Content-Disposition'] = 'attachment; filename="dropkick_labels.csv"'
    return response

def download_h5ad(request):
    file = open('media/dropkick_filter.h5ad', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)

    # decide the file name
    response['Content-Disposition'] = 'attachment; filename="dropkick_filter.h5ad"'
    return response

def download_counts(request):
    file = open('media/dropkick_counts.csv', 'rb')
    response = FileResponse(file)
    
    response['Content-Disposition'] = 'attachment; filename="dropkick_counts.csv"'
    return response

def download_genes(request):
    file = open('media/dropkick_genes.csv', 'rb')
    response = FileResponse(file)
    
    response['Content-Disposition'] = 'attachment; filename="dropkick_genes.csv"'
    return response

def download_sample(request):
    file = open('media/t_4k_small_dropkick_scores.csv', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)
    
    response['Content-Disposition'] = 'attachment; filename="sample_dropkick_scores.csv"'
    return response

def download_qc(request):
    file = open('media/qc_plot.png', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)
    
    response['Content-Disposition'] = 'attachment; filename="qc_plot.png"'
    return response

def download_coef(request):
    file = open('media/coef_plot.png', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)
    
    response['Content-Disposition'] = 'attachment; filename="coef_plot.png"'
    return response

def download_score(request):
    file = open('media/score_plot.png', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)
    
    response['Content-Disposition'] = 'attachment; filename="score_plot.png"'
    return response

def download_all_no_qc(request):
    # Files (local path) to put in the .zip
    # FIXME: Change this (get paths from DB etc)
    filenames = ["media/dropkick_labels.csv", "media/dropkick_counts.csv", "media/dropkick_genes.csv", "media/dropkick_filter.h5ad", "media/coef_plot.png", "media/score_plot.png"]
    
    # Folder name in ZIP archive which contains the above files
    # E.g [thearchive.zip]/somefiles/file2.txt
    # FIXME: Set this to something better
    zip_subdir = "dropkick_output"
    zip_filename = "%s.zip" % zip_subdir
    
    # Open StringIO to grab in-memory ZIP contents
    s = BytesIO()

    # The zip compressor
    zf = zipfile.ZipFile(s, "w")

    for fpath in filenames:
        # Calculate path for file in zip
        fdir, fname = os.path.split(fpath)
        zip_path = os.path.join(zip_subdir, fname)

        # Add file, at correct path
        zf.write(fpath, zip_path)

    # Must close zip for all contents to be written
    zf.close()

    # Grab ZIP file from in-memory, make response with correct MIME-type
    response = HttpResponse(s.getvalue(), content_type = "application/x-zip-compressed")

    # ..and correct content-disposition
    response['Content-Disposition'] = 'attachment; filename=%s' % zip_filename

    return response

def download_all(request):
    # Files (local path) to put in the .zip
    # FIXME: Change this (get paths from DB etc)
    filenames = ["media/dropkick_labels.csv", "media/dropkick_counts.csv", "media/dropkick_genes.csv", "media/dropkick_filter.h5ad", "media/qc_plot.png", "media/coef_plot.png", "media/score_plot.png"]
    
    # Folder name in ZIP archive which contains the above files
    # E.g [thearchive.zip]/somefiles/file2.txt
    # FIXME: Set this to something better
    zip_subdir = "dropkick_output"
    zip_filename = "%s.zip" % zip_subdir
    
    # Open StringIO to grab in-memory ZIP contents
    s = BytesIO()

    # The zip compressor
    zf = zipfile.ZipFile(s, "w")

    for fpath in filenames:
        # Calculate path for file in zip
        fdir, fname = os.path.split(fpath)
        zip_path = os.path.join(zip_subdir, fname)

        # Add file, at correct path
        zf.write(fpath, zip_path)

    # Must close zip for all contents to be written
    zf.close()

    # Grab ZIP file from in-memory, make response with correct MIME-type
    response = HttpResponse(s.getvalue(), content_type = "application/x-zip-compressed")

    # ..and correct content-disposition
    response['Content-Disposition'] = 'attachment; filename=%s' % zip_filename

    return response
