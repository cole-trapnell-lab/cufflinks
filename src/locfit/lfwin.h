#define LFM_EXIT   0
#define LFM_COPY   1
#define LFM_PASTE  2
#define LFM_RUN    3

#define LFM_READA 10
#define LFM_SAVED 11
#define LFM_READD 12
#define LFM_SUMD  13
#define LFM_PLOTD 18

#define LFM_LOCF  20
#define LFM_READF 22
#define LFM_SUMF  23
#define LFM_PRFIT 24

#define LFM_ALPH  70
#define LFM_FIXH  71
#define LFM_APEN  72
#define LFM_DEG0  75
#define LFM_DEG1  76
#define LFM_DEG2  77
#define LFM_DEG3  78

#define LFM_ABOUT 81
#define LFM_INDEX 82
#define LFM_READM 83
#define LFM_WWW   84

#define LFP_ROT   10
#define LFP_STY   11
#define LFP_PS    42
#define LFP_COL   13

#define LFP_XLAB  20
#define LFP_YLAB  21
#define LFP_ZLAB  22
#define LFP_MAIN  23

#define AB_WWW 10

#define CM_LINE 1
#define CM_OK   99

#define RL_ALP  0
#define RL_ALPV 1
#define RL_H    2
#define RL_HV   3
#define RL_PEN  4
#define RL_PENV 5
#define RL_DEG  10
#define RL_FORM 20
#define RL_FAMY 21
#define RL_QUAS 22
#define RL_ROBU 23
#define RL_FIT  98
#define RL_OK   99

#define RP_VS 1
#define RP_HS 2
#define RP_AUT  3
#define RP_DRAW 98
#define RP_OK   99

#define PD_X 1
#define PD_Y 2
#define PD_Z 3
#define PD_DRAW 10
#define PD_ADD  11
#define PD_WIN  12

#define PS_FIL 1
#define PS_DR  8
#define PS_CA  9
#define PS_H 10
#define PS_W 11

#define SC_COL 1
#define SC_SCO 2
#define SC_DR  8
#define SC_OK  9

#define VN_VN 1
#define VN_SA 2
#define VN_RF 98
#define VN_CA 99

#define BP_ALP 1
#define BP_ALV 2
#define BP_AUT 3
#define BP_FIT 4
#define BP_EX 99

#define GR_CM 10
#define GR_ST 11

#define LB_LAB  10
#define LB_DRAW 11

#define LD_QUIT 99

/* about.c */
extern void AboutDlg();

/* devwin.c */
extern void getwinsize(), GetFontInfo();

/* dlgraph.c */
extern void GStyleDlg(), LabelDlg(), PostDlg(), RotateDlg(), SetColDlg();

/* winfile.c */
extern void ReadFileDlg(), ReadDataDlg(), SaveDataDlg(), RunDlg();
extern void ReadFitDlg();

/* windlg.c */
extern void BandDlg(), LocfitDlg(), PlotDataDlg(), wdispatch();
extern int LFDefDlgProc();
