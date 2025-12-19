project(
  'rtiow',
  'cpp',
  meson_version : '>= 1.3.0',
  version : '0.1',
  default_options : ['c_std=c++20', 'warning_level=3'],
)

cc = meson.get_compiler('cpp')
dependencies = [
  cc.find_library('m', required : true)
]

exe = executable(
  'rtiow',
  'rtiow.cpp',
  dependencies : dependencies,
  install : true,
)

test('rtiow', exe)
