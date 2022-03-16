from django.shortcuts import render, get_object_or_404, redirect
from dropkickApp.models import MyFile, CustomParam
from django.views import generic
from .forms import UploadFileForm, CheckboxForm, CustomForm, ScoreForm
from django.http import HttpResponse, StreamingHttpResponse, FileResponse
from django.core.files.storage import FileSystemStorage
import csv
import os
import zipfile
import datetime
from io import BytesIO
from django.core.exceptions import ValidationError
from django.contrib import messages

import scanpy as sc; sc.set_figure_params(color_map="viridis", frameon=False)
import dropkick as dk
import matplotlib.pyplot as plt; plt.switch_backend("Agg")
import io, base64, urllib
import numpy as np
import pandas as pd


def qc_plot(adata, instance):
    # plot QC metrics
    adata = dk.recipe_dropkick(adata, n_hvgs=None, X_final="raw_counts")
    qc_plt = dk.qc_summary(adata)
    
    # display chart
    buf = io.BytesIO()
    qc_plt.savefig(buf, format = 'png')
    qc_plt.savefig('media/qc_plot_' + instance.name + '_' + str(instance.id) + '.png')
    buf.seek(0)
    string = base64.b64encode(buf.read())
    uri = urllib.parse.quote(string)
    return uri

def labels(adata, min_genes, mito_names, n_ambient, n_hvgs, thresh_methods, alphas, max_iter, seed, instance):
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
    coef_plt.savefig('media/coef_plot_' + instance.name + '_' + str(instance.id) + '.png')
    buf_coef.seek(0)
    string_coef = base64.b64encode(buf_coef.read())
    uri_coef = urllib.parse.quote(string_coef)
    
    # display score plot
    adata_score = dk.recipe_dropkick(adata, n_hvgs=None, verbose=False, filter=True, min_genes=50)
    score_plt = dk.score_plot(adata_score)
    buf_score = io.BytesIO()
    score_plt.savefig(buf_score, format = 'png')
    score_plt.savefig('media/score_plot_' + instance.name + '_' + str(instance.id) + '.png')
    buf_score.seek(0)
    string_score = base64.b64encode(buf_score.read())
    uri_score = urllib.parse.quote(string_score)
    
    return uri_score, uri_coef

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
        
def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[:-len(suffix)]
    return input_string

def index(request):
    """View function for home page of site."""
    form = CustomForm(request.POST or None, initial = {'min_genes': 50, 'mito_names': 'mt'})
    model = CustomParam
    # upload file
    if request.method == 'POST':
        if 'document' in request.FILES:
            if form.is_valid():
                instance = form.save()
                request.session['id'] = instance.id
#                 global getID
#                 def getID():
#                     return instance.id

                if instance.qc_plot or instance.dropkick:
                    uploaded_file = request.FILES['document']
                    if uploaded_file.name.endswith('.csv'):
                        fs = FileSystemStorage()
                        fs.save(uploaded_file.name, uploaded_file)
                        instance.name = remove_suffix(uploaded_file.name, '.csv')
                        os.rename('media/' + uploaded_file.name, 'media/' + instance.name + '_' + str(instance.id) + '.csv')
                        instance.csv_bool = True
                        instance.save()
                        param_assignment(instance)
                        return redirect(process)
                    elif uploaded_file.name.endswith('.h5ad'):
                        fs = FileSystemStorage()
                        fs.save(uploaded_file.name, uploaded_file)
                        instance.name = remove_suffix(uploaded_file.name, '.h5ad')
                        os.rename('media/' + uploaded_file.name, 'media/' + instance.name + '_' + str(instance.id) + '.h5ad')
                        instance.h5ad_bool = True
                        instance.save()
                        param_assignment(instance)
                        return redirect(process)
                    elif uploaded_file.name.endswith('.tsv'):
                        fs = FileSystemStorage()
                        fs.save(uploaded_file.name, uploaded_file)
                        instance.name = remove_suffix(uploaded_file.name, '.tsv')
                        os.rename('media/' + uploaded_file.name, 'media/' + instance.name + '_' + str(instance.id) + '.tsv')
                        instance.tsv_bool = True
                        instance.save()
                        param_assignment(instance)
                        return redirect(process)
                    else:
                        messages.error(request,'Please upload a file of CSV, H5AD, or TSV type')
                else:
                    messages.error(request, 'Please select an action to run.')
                
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
    cur_id = request.session['id']
    instance = model.objects.filter(id=cur_id)[0]
    if instance.csv_bool:
        adata = sc.read('media/' + instance.name + '_' + str(instance.id) + '.csv')
    elif instance.h5ad_bool:
        adata = sc.read('media/' + instance.name + '_' + str(instance.id) + '.h5ad')
    elif instance.tsv_bool:
        adata = sc.read('media/' + instance.name + '_' + str(instance.id) + '.tsv')
    
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
        context['qc_plot'] = qc_plot(adata, instance)
        request.session['qc_uri'] = context['qc_plot']
        
    if instance.dropkick:
        # filter checkbox was checked

        # run dropkick
        context['counts_text'] = 'Droplets Inventory'
        context['score_text'] = 'Score Plot'
        context['coef_text'] = 'Coefficient Plot'
        context['labels_text'] = 'Dropkick Labels'
        
        alphas_list = instance.alphas.split(',')
        alphas = [float(x) for x in alphas_list]


        context['score_plot'], context['coef_plot'] = labels(
            adata, instance.min_genes, instance.mito_names, instance.n_ambient, instance.n_hvgs, instance.thresh_methods, alphas,
            instance.max_iter, instance.seed, instance)
        request.session['score_uri'] = context['score_plot']
        request.session['coef_uri'] = context['coef_plot']
        
        context['score_thresh'] = instance.score_thresh

        adata.obs['dropkick_label'] = adata.obs['dropkick_score'] > instance.score_thresh

        context['counts_false'] = adata.obs['dropkick_label'].value_counts()[0]
        context['counts_true'] = adata.obs['dropkick_label'].value_counts()[1]

        # convert dataframe to csv
        adata.obs[["dropkick_score","dropkick_label"]].to_csv('media/dropkick_labels_' + instance.name + '_' + str(instance.id) + '.csv')

        # convert to h5ad file
        adata.write('media/dropkick_filter_' + instance.name + '_' + str(instance.id) + '.h5ad', compression='gzip')

        # output counts and genes matrices
        # output counts csv
        #data_out = adata[df.obs['dropkick_label']==True]
        #data = pd.DataFrame(data_out.X.toarray())
        #data.to_csv('media/dropkick_counts_' + instance.name + '_' + str(instance.id) + '.csv', header=False, index=False)
        
        # output genes csv
        adata.var[["pct_dropout_by_counts","ambient","dropkick_coef"]].to_csv('media/dropkick_genes_' + instance.name + '_' + str(instance.id) + '.csv', header=False, index=False)

    return render(request, 'process.html', context)


def calc_score_thresh(request):
    context = {
        'score_thresh': None, 'title': None, 'qc_text': None, 'counts_text': None, 'counts_false': None, 'counts_true': None,
    }
    form = ScoreForm(request.POST or None)
    model = CustomParam
    cur_id = request.session['id']
    instance = model.objects.filter(id=cur_id)[0]
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
        context['qc_plot'] = request.session['qc_uri']
        
    if instance.dropkick:
        # filter checkbox was checked

        # run dropkick
        context['counts_text'] = 'Droplets Inventory'
        context['score_text'] = 'Score Plot'
        context['score_plot'] = request.session['score_uri']
        context['coef_text'] = 'Coefficient Plot'
        context['coef_plot'] = request.session['coef_uri']
        context['labels_text'] = 'Dropkick Labels'
    
        score_thresh = instance.score_thresh
        context['score_thresh'] = score_thresh
                
        adata = sc.read('media/dropkick_filter_' + instance.name + '_' + str(instance.id) + '.h5ad')
                
        adata.obs['dropkick_label'] = adata.obs['dropkick_score'] > score_thresh

        context['counts_false'] = adata.obs['dropkick_label'].value_counts()[0]
        context['counts_true'] = adata.obs['dropkick_label'].value_counts()[1]
        

        # convert dataframe to csv
        adata.obs[["dropkick_score","dropkick_label"]].to_csv('media/dropkick_labels_' + instance.name + '_' + str(instance.id) + '.csv')

        # convert to h5ad file
        adata.write('media/dropkick_filter_' + instance.name + '_' + str(instance.id) + '.h5ad', compression='gzip')
    return render(request, 'score_thresh.html', context)

def download_csv(request):
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    file = open('media/dropkick_labels_' + instance.name + '_' + str(instance.id) + '.csv', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)

    # decide the file name
    new_filename = 'dropkick_labels_' + instance.name + '_' + str(instance.id);
    response['Content-Disposition'] = 'attachment; filename=%s' % new_filename
    return response

def download_h5ad(request):
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    file = open('media/dropkick_filter_' + instance.name + '_' + str(instance.id) + '.h5ad', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)

    # decide the file name
    new_filename = 'dropkick_filter_' + instance.name + '_' + str(instance.id);
    response['Content-Disposition'] = 'attachment; filename=%s' % new_filename
    return response

def download_counts(request):
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    file = open('media/dropkick_counts_' + instance.name + '_' + str(instance.id) + '.csv', 'rb')
    response = FileResponse(file)
    
    new_filename = 'dropkick_counts_' + instance.name + '_' + str(instance.id) + '.csv';
    response['Content-Disposition'] = 'attachment; filename=%s' % new_filename
    return response

def download_genes(request):
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    file = open('media/dropkick_genes_' + instance.name + '_' + str(instance.id) + '.csv', 'rb')
    response = FileResponse(file)
    
    new_filename = 'dropkick_genes_' + instance.name + '_' + str(instance.id);
    response['Content-Disposition'] = 'attachment; filename=%s' % new_filename
    return response

def download_sample(request):
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    file = open('dropkickApp/static/t_4k_small_dropkick_scores.csv', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)
    
    response['Content-Disposition'] = 'attachment; filename="sample_dropkick_scores.csv"'
    return response

def download_qc(request):
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    file = open('media/qc_plot_' + instance.name + '_' + str(instance.id) + '.png', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)
    
    new_filename = 'qc_plot_' + instance.name + '_' + str(instance.id);
    response['Content-Disposition'] = 'attachment; filename=%s' % new_filename
    return response

def download_coef(request):
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    file = open('media/coef_plot_' + instance.name + '_' + str(instance.id) + '.png', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)
    
    new_filename = 'coef_plot_' + instance.name + '_' + str(instance.id);
    response['Content-Disposition'] = 'attachment; filename=%s' % new_filename
    return response

def download_score(request):
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    file = open('media/score_plot_' + instance.name + '_' + str(instance.id) + '.png', 'rb') # Read the file in binary mode, this file must exist
    response = FileResponse(file)
    
    new_filename = 'score_plot_' + instance.name + '_' + str(instance.id);
    response['Content-Disposition'] = 'attachment; filename=%s' % new_filename
    return response

def download_all_no_qc(request):
    # Files (local path) to put in the .zip
    # FIXME: Change this (get paths from DB etc)
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    filenames = ['media/dropkick_labels_' + instance.name + '_' + str(instance.id) + '.csv', 'media/dropkick_genes_' + instance.name + '_' + str(instance.id) + '.csv', 'media/dropkick_filter_' + instance.name + '_' + str(instance.id) + '.h5ad', 'media/coef_plot_' + instance.name + '_' + str(instance.id) + '.png', 'media/score_plot_' + instance.name + '_' + str(instance.id) + '.png']
    
    # Folder name in ZIP archive which contains the above files
    # E.g [thearchive.zip]/somefiles/file2.txt
    # FIXME: Set this to something better
    zip_subdir = "dropkick_output_" + instance.name
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
    cur_id = request.session['id']
    instance = CustomParam.objects.filter(id=cur_id)[0]
    filenames = ['media/dropkick_labels_' + instance.name + '_' + str(instance.id) + '.csv', 'media/dropkick_genes_' + instance.name + '_' + str(instance.id) + '.csv', 'media/dropkick_filter_' + instance.name + '_' + str(instance.id) + '.h5ad', 'media/qc_plot_' + instance.name + '_' + str(instance.id) + '.png', 'media/coef_plot_' + instance.name + '_' + str(instance.id) + '.png', 'media/score_plot_' + instance.name + '_' + str(instance.id) + '.png']
    
    # Folder name in ZIP archive which contains the above files
    # E.g [thearchive.zip]/somefiles/file2.txt
    # FIXME: Set this to something better
    zip_subdir = "dropkick_output_" + instance.name
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
