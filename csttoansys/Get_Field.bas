Option Explicit

'#include "vba_globals_all.lib"

Public Type coordinates
  x As Double
  y As Double
  z As Double
End Type

'Fixed Length of Array, eventually make it later dynamic
Const MAX_NP = 10000000


Sub Main

  Dim CST_nNode(MAX_NP) As Long
  Dim CST_points(MAX_NP) As coordinates
  Dim count As Long, ip As Long
  Dim surface_x As Double, surface_y As Double, surface_z As Double
  Dim res1 As Object
  Set res1 = Result3D("^field_se_phi.m3d")
  Dim info As String
  info = "Choose the output file. Default is PRNSOL.lis."
  Dim outdir As String
  Begin Dialog UserDialog 600,112,"Extract potentials into the given directory as PRNSOL.lis",.DialogFunc ' %GRID:10,7,1,1
    Text 20,21,90,14,"Output File",.Text1
    TextBox 100,21,390,21,.Outfile
    PushButton 500,21,90,21,"Browse",.Browseinputdir
    OKButton 410,84,90,21
    CancelButton 500,84,90,21
    Text 20,49,500,21,info,.Text3
  End Dialog
  Dim dlg As UserDialog
  If (Dialog(dlg) = 0) Then Exit All

  Dim outfile As String
  outfile = dlg.Outfile

  Open  outfile For Output As #2

  '---Loop for mesh points
  With Mesh
  Plot3DPlotsOn2DPlane False
  VectorPlot3D.Reset
  For ip = 0 To .GetNp-1
    Print #2, "  "  + CStr(ip+1) + "  "  +  res1.GetXRe(ip)
    count = count+1
  Next ip
  End With
  Close #2

End Sub
'-------------------------------------------------------------------------------------
Function dialogfunc(DlgItem$, Action%, SuppValue%) As Boolean

  Dim Extension As String, projectdir As String, filename As String
    Select Case Action%
    Case 1 ' Dialog box initialization
    Case 2 ' Value changing or button pressed
        Select Case DlgItem
          Case "Browseinputdir"
'                projectdir = GetProjectPath("Root")
                projectdir = "N:\4all\xxl\zenker\flc_xxl\CST\Export2Garfield\"
                filename = GetFilePath("PRNSOL.lis",, projectdir, "Click On Any file In the folder where the PRNSOL.lis should be created", 2)
                If (filename <> "") Then
                    DlgText "Outfile", filename
                End If
          dialogfunc = True
        End Select
    Case 3 ' TextBox or ComboBox text changed
    Case 4 ' Focus changed
    Case 5 ' Idle
    Case 6 ' Function key
    End Select
End Function

