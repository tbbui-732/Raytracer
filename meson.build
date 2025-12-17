project(
  'RTIOW',
  'c',
  meson_version : '>= 1.3.0',
  version : '0.1',
  default_options : ['c_std=c99', 'warning_level=3'],
)

dependencies = [
]

exe = executable(
  'RTIOW',
  'rtiow.c',
  dependencies : dependencies,
  install : true,
)

test('rtiow', exe)
