# make Apple Silicon app
#
# rm -rf build dist
# python setup_m1.py py2app
# mkdir -p dist/dmg
cp -r "dist/demo.app" dist/dmg
# If the DMG already exists, delete it.
test -f "dist/demo.dmg" && rm "dist/demo.dmg"
# get create-dmg from homebrew
create-dmg \
--volname "demo" \
--volicon "icon.icns" \
--window-pos 200 120 \
--window-size 600 300 \
--icon-size 64 \
--icon "demo.app" 128 128 \
--hide-extension "demo.app" \
--app-drop-link 425 120 \
"dist/demo.dmg" \
"dist/dmg/"

#hdiutil create -srcFolder dist -o dmg/neurodemo

