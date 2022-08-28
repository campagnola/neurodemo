# make Apple Silicon app
#
# The next 3 lines can be skipped if you
# do not need to update the code
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
--window-pos 300 300 \
--window-size 400 300 \
--icon-size 32 \
--icon "demo.app" 100 100 \
--hide-extension "demo.app" \
--app-drop-link 250 100 \
"dist/demo_M1.dmg" \
"dist/dmg/"

#hdiutil create -srcFolder dist -o dmg/neurodemo

