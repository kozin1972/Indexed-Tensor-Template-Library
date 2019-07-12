from distutils.core import setup, Extension

iTTL_path='/d/workspace/tensor_11'
l1_procs_module_path=iTTL_path+'/l1_procs_module.cpp'
BLAS_library='openblas'

#print l1_procs_module_path

module1 = Extension('l1_procs',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = [iTTL_path],
#                    include_dirs = ['/d/workspace/tensor_11'],
                    libraries = [BLAS_library],
#                    library_dirs = ['/usr/local/lib'],
                    sources = [l1_procs_module_path])
#                    sources = ['/d/workspace/tensor_11/l1_procs_module.cpp'])

setup (name = 'iTTL_procs',
       version = '1.0',
       description = 'Indexed Tensor Template Library subroutines demo package',
       author = 'Alexey Kozin',
#       author_email = '@',
       url = 'https://github.com/kozin1972/Indexed-Tensor-Template-Library',
       long_description = '''
This is really just a demo package.
''',
       ext_modules = [module1])