# -*- mode: python ; coding: utf-8 -*-


block_cipher = None

datafiles = [("neurodemo_venv/Lib/site-packages/pyqtgraph/graphicsItems/ViewBox/axisCtrlTemplate_pyqt6.py", 'pyqtgraph'),
    ]
a = Analysis(
    ['setup_win.py'],
    pathex=[],
    binaries=[],
    datas=datafiles,
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='setup_win',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    windowed=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_name='NeuroDemo_1_exe'
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
