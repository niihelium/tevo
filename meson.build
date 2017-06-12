project('tevo', 'fortran')

#
# Compiler flags
#
dialect = ''
bounds =  ''
if meson.get_compiler('fortran').get_id() == 'gcc'
    dialect = [ '-ffree-form','-fimplicit-none']
    bounds = '-fbounds-check'
endif
if meson.get_compiler('fortran').get_id() == 'intel'
    dialect = ['-stand f08','-free','-implicitnone']
    bounds = '-check bounds'
endif
if meson.get_compiler('fortran').get_id() == 'pgi'
    dialect = [ '-Mfreeform','-Mdclchk',
                '-Mstandard','-Mallocatable=03']
    bounds = '-C'
endif

add_global_arguments(dialect, language : 'fortran')
add_global_arguments(bounds, language : 'fortran')

sources = ['src/main.f08',
           'src/reader.f08',
           'src/spectre.f08',
           'src/calculation.f08',
           'src/constants.f08',
           'src/dvode_f90.f90',
           'src/cooler.f08']

executable('tevo', sources)
