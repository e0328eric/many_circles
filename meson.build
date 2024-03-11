project(
  'many_circles',
  'c',
  license: ['MIT'],
  version: '0.1.0',
  meson_version: '>=1.0.0',
  default_options: [
      'buildtype=release',
      'c_std=gnu11',
      'warning_level=1',
  ]
)

gsl = dependency('gsl')

sources = files(
    'src/main.c',
    'src/extended_bessel.c',
    'src/circle.c',
    'src/single_layer.c',
    'src/single_circle_solv.c',
)

executable('many_circles', sources,
  dependencies: [gsl],
)
