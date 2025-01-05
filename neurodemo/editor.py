import sys
from . import qt


class Editor(qt.QPlainTextEdit):
    def __init__(self, parent=None, language=None):
        qt.QPlainTextEdit.__init__(self, parent)

    def setText(self, text):
        self.setPlainText(text)

    def text(self):
        return str(self.toPlainText()).encode('UTF-8')

    def __getattr__(self, name):
        return lambda: None



if __name__ == "__main__":
    app = qt.QApplication(sys.argv)
    editor = Editor()
    editor.show()
    editor.setText(open(sys.argv[0]).read())
    editor.resize(800, 800)
    app.exec_()
