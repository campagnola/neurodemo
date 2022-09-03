# If the DMG already exists, delete it.
test -f "dist/demo_x86_64.dmg" && rm "dist/demo_x86_64.dmg"
sleep 3
test -f "dist/demo.dmg" && rm "dist/demo.dmg"

# make x86_64 app
#
# rm -rf build dist
# python setup_x86_64.py py2app
# mkdir -p dist/dmg
sleep 3

cp -r "dist/demo.app" dist/dmg

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
--no-internet-enable \
"dist/demo_x86_64.dmg" \
"dist/dmg/"

#hdiutil create -srcFolder dist -o dmg/neurodemo

