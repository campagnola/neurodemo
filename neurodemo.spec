# -*- mode: python ; coding: utf-8 -*-

added_files = [
    ( './neurodemo/images/channel.svg', 'images' ),
    ( './neurodemo/images/channel2.svg', 'images' ),
    ( './neurodemo/images/cell.svg', 'images' ),
    ( './neurodemo/images/pipette.svg', 'images' )
]

a = Analysis(
    ['neurodemo.py'],
    pathex=[],
    binaries=[],
    datas=added_files,
    hiddenimports=[ 'PyQt6.QtSvgWidgets' ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)
splash = Splash(
    'splash_screen.png',
    binaries=a.binaries,
    datas=a.datas,
    text_pos=None,
    text_size=12,
    minify_script=True,
    always_on_top=False,
)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    splash,
    splash.binaries,
    [],
    name='neurodemo',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
