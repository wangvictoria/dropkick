from django import forms
from .models import MyFile, CustomParam, Contact
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from django.core.validators import DecimalValidator, validate_integer
from django.utils.datastructures import MultiValueDict

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

class UploadFileForm(forms.ModelForm):
    class Meta:
        model = MyFile
        fields = ('name', 'upload',)

class CheckboxForm(forms.Form):
#     def validate_int(value):
#         if not(value.isdigit()):
#             raise ValidationError(
#                 _('%(value)s is not an integer'),
#                 params={'value': value},
#             )
    def validate_integer(value):
        if not(value.isnumeric()):
            raise ValidationError(_('%(value)s is not a non-negative integer.'),
            params={'value': value},)
        else:
            return value
            
    def IntegerValidator(value):
        if not(validate_integer(value)):
            messages.error(value,'Please enter a non-negative integer.')
    def validate_dec(value):
        try:
            float(element)
        except ValueError:
            raise ValidationError(
                _('%(value)s is not a float'),
                params={'value': value},
            )

    
    # qc plot and/or filter
    qc_plot = forms.BooleanField(label="qc_plot", required=False)
    dropkick = forms.BooleanField(label="dropkick", required=False)
    
    # default or custom settings
    default = forms.BooleanField(label="default", required=False)
    custom = forms.BooleanField(label="custom", required=False)
    
    # parameters
    min_genes = forms.CharField(max_length=10, label="min_genes", required=False, empty_value='50', initial='50')
    mito_names = forms.CharField(max_length=100, label="mito_names", required=False, empty_value='^mt-|^MT-', initial='^mt-')
    n_ambient = forms.CharField(max_length=10, label="n_ambient", required=False, empty_value='10', initial='10')
    n_hvgs = forms.CharField(max_length=10, label="n_hvgs", required=False, empty_value='2000', initial='2000')
    #metrics = forms.CharField(max_length=100, label="metrics", required=False, empty_value='arcsinh_n_genes_by_counts', initial='arcsinh_n_genes_by_counts')
    thresh_methods = forms.ChoiceField(label='thresh_methods', choices=THRESH_METHODS, initial='multiotsu')
    score_thresh = forms.CharField(label='score_thresh', required=False, empty_value='0.5', initial='0.5')
    #directions = forms.ChoiceField(label="directions", choices=DIRECTIONS, initial='above')
    alphas = forms.CharField(max_length=10, label="alphas", required=False, empty_value='0.1', initial='0.1') # FIX FOR MULTIPLE ALPHAS
    max_iter = forms.CharField(max_length=10, label="max_iter", required=False, empty_value='2000', initial='2000')
    #n_jobs = forms.CharField(max_length=10, label="n_jobs", required=False, empty_value='2', initial='2')
    seed = forms.CharField(max_length=10, label="seed", required=False, empty_value='18', initial='18')
    
    def cleaned_data(self):
        cd = self.cleaned_data
 
        if not(validate_integer(cd.get('min_genes'))):
            self._errors['min_genes'] = self.error_class([
                'Minimum 5 characters required'])
        if validate_integer(cd.get('n_ambient', None)):
            self.add_error('n_ambient', 'Please enter a non-negative integer.')
        if validate_integer(cd.get('n_hvgs, None')):
            self.add_error('n_hvgs', 'Please enter a non-negative integer.')
        # HOW TO VALIDATE ALPHAS?
        validate_integer(cd.get('max_iter', None))
        validate_integer(cd.get('seed', None))
 
        # return any errors if found
        return cd
    
    
#     def clean(self):
#         # data from the form is fetched using super function
#         super(CheckboxForm, self).clean()
         
#         # extract the username and text field from the data
#         username = self.cleaned_data.get('username')
#         text = self.cleaned_data.get('text')
 
#         # conditions to be met for the username length
#         if len(username) < 5:
#             self._errors['username'] = self.error_class([
#                 'Minimum 5 characters required'])
#         if len(text) <10:
#             self._errors['text'] = self.error_class([
#                 'Post Should Contain a minimum of 10 characters'])
 
#         # return any errors if found
#         return self.cleaned_data

class CustomForm(forms.ModelForm):
#     def clean(self):
#         min_genes = self.cleaned_data['min_genes']
#         if not min_genes:
#             min_genes = 50

#         return min_genes
    
    class Meta:
        model = CustomParam
        fields = ['qc_plot', 'dropkick', 'default', 'custom', 'min_genes', 'mito_names', 'n_ambient', 'n_hvgs', 'thresh_methods', 'score_thresh', 'alphas', 'max_iter', 'seed']
        default={'min_genes': 50}
        #labels = {'qc_plot': _('QC Plot'), 'min_genes': _('Minimum genes')}
        #help_texts = {'min_genes': _('Enter a non-negative integer.')}

class ScoreForm(forms.ModelForm):
    class Meta:
        model = CustomParam
        fields = ['score_thresh',]
        
class ContactForm(forms.ModelForm):
    class Meta:
        model = Contact
        fields = [
            'name',
            'email',
            'comments',
        ]