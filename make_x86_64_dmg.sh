# If the DMG already exists, delete it.
test -f "dist/demo_x86_64.dmg" && rm "dist/demo_x86_64.dmg"
sleep 3
test -f "dist/neurodemo.dmg" && rm "dist/neurodemo.dmg"

# make x86_64 app
#
# rm -rf build dist
# python setup_x86_64.py py2app
# mkdir -p dist/dmg
sleep 3

cp -r "dist/neurodemo.app" dist/dmg

# get create-dmg from homebrew
create-dmg \
--volname "neurodemo" \
--volicon "icon.icns" \
--window-pos 300 300 \
--window-size 400 300 \
--icon-size 32 \
--icon "neurodemo.app" 100 100 \
--hide-extension "neurodemo.app" \
--app-drop-link 250 100 \
--no-internet-enable \
"dist/demo_x86_64.dmg" \
"dist/dmg/"

# hdiutil create -srcFolder dist -o dmg/neurodemo

