import ROOT

def makeText(text, x, y, font=43, size=12, align=12, color=1):
    latex = ROOT.TLatex(3.5, 24, text)
    latex.SetNDC()
    latex.SetTextAlign(align)
    latex.SetTextFont(font)
    latex.SetTextSize(size)
    latex.SetTextColor(color)
    latex.SetX(x)
    latex.SetY(y)
    return latex

def makeLogo(text, x, y):
    latex = ROOT.TLatex(3.5, 24, text)
    latex.SetNDC()
    latex.SetTextFont(43)
    latex.SetTextSize(12)
    latex.SetTextColor(1)
    return latex

def makeSubtitle(text, x, y):
    latex = ROOT.TLatex(3.5, 24, text)
    latex.SetNDC()
    latex.SetTextFont(43)
    latex.SetTextSize(12)
    latex.SetTextColor(1)
    return latex
