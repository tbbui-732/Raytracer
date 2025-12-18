project(
  'rtiow',
  'c',
  meson_version : '>= 1.3.0',
  version : '0.1',
  default_options : ['c_std=c99', 'warning_level=3'],
)

cc = meson.get_compiler('c')
dependencies = [
  cc.find_library('m', required : true)
]

exe = executable(
  'rtiow',
  'rtiow.c',
  dependencies : dependencies,
  install : true,
)

test('rtiow', exe)
