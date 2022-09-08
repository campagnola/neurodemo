# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files, collect_all

block_cipher = None
datafiles = []
binaries = []
datas = []

a = Analysis(
    ['demo.py'],
    pathex=[],
    binaries=binaries,
    datas = [('neurodemo\\images', 'images'),
        ],
    hiddenimports=[
        'pyqtgraph.graphicsItems.ViewBox.axisCtrlTemplate_pyqt6',
        'pyqtgraph.graphicsItems.PlotItem.plotConfigTemplate_pyqt6',
        'pyqtgraph.imageview.ImageViewTemplate_pyqt6',
        'pyqtgraph.console.template_pyqt6',
        'ctypes',
        'numpy',
        'scipy',
        # 'c:\api-ms-win-shcore-scaling-l1-1-1.dll',
        ],
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
    name='NeuroDemo',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    windowed=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_name='NeuroDemo.exe',
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
