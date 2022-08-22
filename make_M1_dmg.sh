# make Apple Silicon app
#
# The next 3 lines can be skipped if you
# do not need to update the code
rm -rf build dist
python setup_m1.py py2app
mkdir -p dist/dmg

cp -r "dist/demo_M1.app" dist/dmg
# If the DMG already exists, delete it.
test -f "dist/demo_M1.dmg" && rm "dist/demo_M1.dmg"
# get create-dmg from homebrew
create-dmg \
--volname "demo" \
--volicon "icon.icns" \
--window-pos 200 120 \
--window-size 600 300 \
--icon-size 64 \
--icon "demo.app" 128 128 \
--hide-extension "demo_M1.app" \
--app-drop-link 425 120 \
"dist/demo_M1.dmg" \
"dist/dmg/"

#hdiutil create -srcFolder dist -o dmg/neurodemo

