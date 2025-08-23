## @package EnableTesterGeotechnical
# The EnableTesterGeotechnical module can be used in all Geotechnical materials by simply adding 1 import line to
# enable the tester
#
# from opensees.physical_properties.utils.tester.EnableTesterGeotechnical import *

from opensees.physical_properties.utils.tester.TesterGeotechnical import GeotechnicalTraits, TesterGeotechnicalWidget

class TesterGeotechnicalGuiGlobals:
    # stores a reference to the gui generated for this object
    gui = None

def __removeGui():
    if TesterGeotechnicalGuiGlobals.gui is not None:
        TesterGeotechnicalGuiGlobals.gui.setParent(None)
        TesterGeotechnicalGuiGlobals.gui.deleteLater()
        TesterGeotechnicalGuiGlobals.gui = None

def onEditorClosing(editor, xobj):
    __removeGui()

def onEditFinished(editor, xobj):
    if TesterGeotechnicalGuiGlobals.gui is not None:
        TesterGeotechnicalGuiGlobals.gui.onEditFinished()

def onEditBegin(editor, xobj):
    __removeGui()
    TesterGeotechnicalGuiGlobals.gui = TesterGeotechnicalWidget(GeotechnicalTraits.D3, editor, xobj)
