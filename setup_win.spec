# -*- mode: python ; coding: utf-8 -*-


block_cipher = None

datafiles = []
#("neurodemo_venv\Lib\site-packages\pyqtgraph\graphicsItems\ViewBox\\axisCtrlTemplate_pyqt6.py", 'pyqtgraph'),
#            ("neurodemo_venv\Lib\site-packages\pyqtgraph\graphicsItems\PlotItem\\plotConfigTemplate_pyqt6.py", 'pyqtgraph'),
#            ("neurodemo_venv\Lib\site-packages\pyqtgraph\imageview\\ImageViewTemplate_pyqt6.py",'pyqtgraph'),
#            ("neurodemo_venv\Lib\site-packages\pyqtgraph\console\template_pyqt6.py", "pyqtgraph"),
#            ]
a = Analysis(
    ['demo.py'],
    pathex=[],
    binaries=[],
    datas=datafiles,
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
