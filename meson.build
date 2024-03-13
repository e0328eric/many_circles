project(
  'many_circles',
  'c',
  license: ['MIT'],
  version: '0.1.0',
  meson_version: '>=1.0.0',
  default_options: [
      'buildtype=release',
      'c_std=gnu11',
      'warning_level=2',
  ]
)

gsl = dependency('gsl')
complex_bessel = declare_dependency(
  link_args: [
    '-lcomplex_bessel',
    '-L/home/almagest/Github/Many_Circles/complex_bessel/build/lib',
  ],
)

sources = files(
    'src/main.c',
    'src/bessel_func.c',
    'src/circle.c',
    'src/single_layer.c',
    'src/single_circle_solv.c',
)

executable('many_circles', sources,
  dependencies: [gsl, complex_bessel],
)
