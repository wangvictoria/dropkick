import datetime
from django.db import models
from django.http import HttpResponseRedirect
from django.shortcuts import render
#from .forms import ModelFormWithFileField

# Create your models here.

THRESH_METHODS = [
    ('multiotsu', 'multiotsu'),
    ('otsu', 'otsu'),
    ('li', 'li'),
    ('mean', 'mean'),
]

DIRECTIONS = [
    ('above', 'above'),
    ('below', 'below'),
]

def user_directory_path(instance, filename):
    # file will be uploaded to MEDIA_ROOT/user_<id>/<filename>
    return 'user_{0}/{1}'.format(instance.user.id, filename)

class MyFile(models.Model):
    """Model representing a file that the user uploads."""
    name = models.CharField('file name', max_length=500, help_text='Enter file output name.', blank=True)
    upload = models.FileField(upload_to='uploads/', help_text='Only H5AD files are supported.')
    uploaded_at = models.DateTimeField('upload date', auto_now_add=True)
    
    class Meta:
        verbose_name = "file"
        ordering = ['uploaded_at', 'name']
    
    def __str__(self):
        """String for representing the Model object."""
        return self.name
    
    def get_absolute_url(self):
        """Returns the url to access a detail record for this file."""
        return reverse('file-detail', args=[str(self.id)])

class CustomParam(models.Model):
    # unique datetime id
    datetime = models.DateTimeField(unique=True, auto_now_add=True)
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=50, default='sample')
    
    # qc plot and/or filter
    qc_plot = models.BooleanField(verbose_name="qc_plot")
    dropkick = models.BooleanField(verbose_name="dropkick")
    
    # default or custom settings
    default = models.BooleanField(verbose_name="default")
    custom = models.BooleanField(verbose_name="custom")
    
    # type of file
    csv_bool = models.BooleanField(default=False)
    h5ad_bool = models.BooleanField(default=False)
    tsv_bool = models.BooleanField(default=False)
    
    # parameters
    min_genes = models.IntegerField(verbose_name="min_genes", default='50', blank=True, null=True)
    mito_names = models.CharField(max_length=100, verbose_name="mito_names", default='^mt-|^MT-', blank=True, null=True)
    n_ambient = models.IntegerField(verbose_name="n_ambient", default='10', blank=True, null=True)
    n_hvgs = models.IntegerField(verbose_name="n_hvgs", default='2000', blank=True, null=True)
    metrics = models.CharField(max_length=100, verbose_name="metrics", default='arcsinh_n_genes_by_counts', blank=True, null=True)
    thresh_methods = models.CharField(max_length=10, verbose_name="thresh_methods", choices=THRESH_METHODS, default='multiotsu', blank=True, null=True)
    score_thresh = models.DecimalField(decimal_places=5, max_digits=10, verbose_name="score_thresh", default=0.5, blank=True, null=True)
    directions = models.CharField(max_length=10, verbose_name="directions", choices=DIRECTIONS, default='above', blank=True, null=True)
    alphas = models.CharField(max_length=10, verbose_name="alphas", default='0.1', blank=True, null=True)
    max_iter = models.IntegerField(verbose_name="max_iter", default='2000', blank=True, null=True)
    n_jobs = models.IntegerField(verbose_name="n_jobs", default='5', blank=True, null=True)
    seed = models.IntegerField(verbose_name="seed", default='18', blank=True, null=True)
    
    score_thresh = models.DecimalField(max_digits=9,
                              decimal_places=5, default=0.5, blank=True, null=True,
                              error_messages={'required': "Please enter a value between 0 and 1."})
    