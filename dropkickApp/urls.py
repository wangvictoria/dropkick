from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('rescore/', views.calc_score_thresh, name='calc_score_thresh'),
    path('download_csv/', views.download_csv, name='download_csv'),
    path('download_h5ad/', views.download_h5ad, name='download_h5ad'),
    path('download_counts/', views.download_counts, name='download_counts'),
    path('download_genes/', views.download_genes, name='download_genes'),
    path('download_sample/', views.download_sample, name='download_sample'),
    path('download_qc/', views.download_qc, name='download_qc'),
    path('download_coef/', views.download_coef, name='download_coef'),
    path('download_score/', views.download_score, name='download_score'),
    path('download_all_no_qc/', views.download_all_no_qc, name='download_all_no_qc'),
    path('download_all/', views.download_all, name='download_all'),
    #path('upload/', views.upload_file, name='upload'),
    path('process/', views.process, name='process'),
]