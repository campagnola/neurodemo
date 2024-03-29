# -*- mode: python ; coding: utf-8 -*-


block_cipher = None

datafiles =[
#            ("c:\neurodemo\neurodemo_venv\Lib\site-packages\pyqtgraph\graphicsItems\ViewBox\axisCtrlTemplate_pyqt6.py", 'pyqtgraph'),
#            ("c:\neurodemo\neurodemo_venv\Lib\site-packages\pyqtgraph\graphicsItems\PlotItem\plotConfigTemplate_pyqt6.py", 'pyqtgraph'),
#            ("c:\neurodemo\neurodemo_venv\Lib\site-packages\pyqtgraph\imageview\ImageViewTemplate_pyqt6.py",'pyqtgraph'),

            ("neurodemo_venv\Lib\site-packages\pyqtgraph\graphicsItems\ViewBox\\axisCtrlTemplate_pyqt6.py", 'pyqtgraph'),
            ("neurodemo_venv\Lib\site-packages\pyqtgraph\graphicsItems\PlotItem\\plotConfigTemplate_pyqt6.py", 'pyqtgraph'),
            ("neurodemo_venv\Lib\site-packages\pyqtgraph\imageview\\ImageViewTemplate_pyqt6.py",'pyqtgraph'),
            ]

a = Analysis(
    ['test_installer.py'],
    pathex=[],
    binaries=[],
    datas=datafiles, 
    hiddenimports=[
                'pyqtgraph.graphicsItems.ViewBox.axisCtrlTemplate_pyqt6',
                'pyqtgraph.graphicsItems.PlotItem.plotConfigTemplate_pyqt6',
                'pyqtgraph.imageview.ImageViewTemplate_pyqt6',
                'numpy', 'scipy',
             ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["tkinter"],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='test_installer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='test_installer',
)
