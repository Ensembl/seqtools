/*
 *  utilities.h
 *  acedb
 *
 *  Created by Gemma Barson on 05/01/2010.
 *
 */
 
#ifndef _utilities_h_included_
#define _utilities_h_included_

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>

#define UNSET_INT  -1


#define max(a,b)        (((a) > (b)) ? (a) : (b))
#define min(a,b)        (((a) < (b)) ? (a) : (b))


/* Debug logging macros. #define DEBUG to enable debug output. */
#ifdef DEBUG
#define DEBUG_OUT(format, args...) debugLogLevel(0); printf(format, ##args);
#else
#define DEBUG_OUT(format, args...)
#endif

#ifdef DEBUG
#define DEBUG_ENTER(format, args...) debugLogLevel(1); printf("-->"); printf(format, ##args); printf("\n");
#else
#define DEBUG_ENTER(format, args...)
#endif

#ifdef DEBUG
#define DEBUG_EXIT(format, args...) debugLogLevel(-1); printf("<--"); printf(format, ##args); printf("\n");
#else
#define DEBUG_EXIT(format, args...)
#endif



/* Color strings that can be passed to parseBlxColor to create a GdkColor */
#define BLX_YELLOW	      "#ffff00" 
#define BLX_DARK_YELLOW	      "#c0c000"
#define BLX_LIGHT_YELLOW      "#fafad2"
#define BLX_CYAN	      "#00ffff"
#define BLX_DARK_CYAN	      "#008080"
#define BLX_LIGHT_CYAN	      "#80FFF2"
#define BLX_VERY_LIGHT_BLUE   "#c1ccdb"
#define BLX_BLUE	      "#0000ff"
#define BLX_DARK_BLUE	      "#000080"
#define BLX_PALE_STEEL_BLUE   "#78b4f0"
#define BLX_LIGHT_STEEL_BLUE  "#b0c4ff"
#define BLX_STEEL_BLUE	      "#5e76bb"
#define BLX_ROYAL_BLUE	      "#4169e1"
#define BLX_LIGHT_SKY_BLUE    "#87cefa"
#define BLX_SKY_BLUE	      "#a0b8c8"
#define BLX_GREY	      "#bebebe"
#define BLX_DARK_GREY	      "#929292"
#define BLX_LIGHT_GREY	      "#d3d3d3"
#define BLX_BLACK	      "#000000"
#define BLX_WHITE	      "#ffffff"
#define BLX_RED		      "#ff0000"
#define BLX_LIGHTER_RED	      "#ff5353"
#define BLX_LIGHT_RED	      "#ff7373"
#define BLX_DARK_RED	      "#800000"
#define BLX_VERY_DARK_RED     "#A00000"
#define BLX_ORANGE	      "#ffa500"
#define BLX_LIGHT_ORANGE      "#FFCF55"
#define BLX_DUSK_ORANGE	      "#D7BA94"
#define BLX_GREEN	      "#00ff00"
#define BLX_DARK_GREEN	      "#00bb00"
#define BLX_VERY_DARK_GREEN   "#015800"
#define BLX_TURQUOISE	      "#40e0d0"
#define BLX_DARK_TURQUOISE    "#32aca0"
#define BLX_PURPLE	      "#9370db"
#define BLX_PLUM	      "#dda0dd"
#define BLX_THISTLE	      "#d8bfd8"
#define BLX_ORANGE_RED	      "#ff4500"
#define BLX_PALE_GREEN	      "#C1FFC1" 
#define BLX_LAWN_GREEN	      "#7cfc00"
#define BLX_LIGHT_CORAL	      "#f08080"
#define BLX_LIGHT_SALMON      "#ffa07a"
#define BLX_BURLYWOOD	      "#deb887"
#define BLX_TAN		      "#d2b48c"

/* The following are used to define default colors for certain types of features in Blixem.
 * One of several different actual colors from the BlxColor struct may be used depending 
 * on state, e.g. we use a different color if "print colors" (i.e. black and 
 * white mode) is on. */

typedef enum 
  {
    BLXCOL_MIN,		  /* dummy value so that we don't get a zero ID */
  
    BLXCOL_BACKGROUND,	  /* background color of the widgets */
    BLXCOL_REF_SEQ,	  /* default background color for the reference sequence */  
    BLXCOL_MATCH,	  /* background color for an exact match */
    BLXCOL_CONS,	  /* background color for a conserved match */
    BLXCOL_MISMATCH,	  /* background color for a mismatch */
    BLXCOL_EXON_CDS,	  /* background color for an exon (coding region) */
    BLXCOL_EXON_UTR,	  /* background color for an exon (non-coding/untranslated region) */
    BLXCOL_INSERTION,	  /* color for an insertion marker */
    BLXCOL_EXON_START,	  /* color for the start boundary line of an exon */
    BLXCOL_EXON_END,  	  /* color for the end boundary line of an exon */
    BLXCOL_CODON,	  /* color in which to highlight the nucleotides for the currently-selected codon */
    BLXCOL_MET,		  /* background color for MET codons in the three frame translation */
    BLXCOL_STOP,	  /* background color for STOP codons in the three frame translation */
    BLXCOL_GRID_LINE,	  /* color of the gridlines in the big picture grids */
    BLXCOL_GRID_TEXT,	  /* color of the text in the big picture grids */
    BLXCOL_HIGHLIGHT_BOX, /* color of the highlight box in the big picture */
    BLXCOL_PREVIEW_BOX,	  /* color of the preview box in the big picture */
    BLXCOL_MSP_LINE,	  /* color of the MSP lines in the big picture */
    BLXCOL_SNP,		  /* background color for SNPs */
    BLXCOL_GROUP,	  /* default highlight color for generic groups */
    BLXCOL_MATCH_SET,	  /* default highlight color for the special match-set group */
    BLXCOL_EXON_FILL_CDS, /* fill color for an exon in the big picture (coding region) */
    BLXCOL_EXON_FILL_UTR, /* fill color for an exon in the big picture (non-coding/untranslated region) */
    BLXCOL_EXON_LINE_CDS, /* line color for an exon in the big picture (coding region) */
    BLXCOL_EXON_LINE_UTR, /* line color for an exon in the big picture (non-coding/untranslated region) */
    BLXCOL_UNALIGNED_SEQ, /* color in which to show additional sequence in the match that is not part of the alignment */

    BLXCOL_NUM_COLORS
  } BlxColorId;

typedef struct _BlxColor
  {
    BlxColorId id;		  /* unique ID for this color */
    char *name;			  /* meaningful name for the color e.g. "Match" */
    char *desc;			  /* meaningful description for what the color is used for e.g. "background color for exact matches" */
    gboolean transparent;	  /* if this is true, the colors below are not specified and the background color should be used instead */
    
    GdkColor normal;		  /* the color in normal operation */
    GdkColor selected;		  /* the color in a selected state */
    GdkColor print;		  /* the color used for printing */
    GdkColor printSelected;	  /* the selected-state color used for printing */
  } BlxColor;



/* This struct holds generic callback information. It can be stored as a widget
 * property and called on the widget on request (e.g. by a dialog when it applies changes). */
typedef struct _CallbackData
  {
    GtkCallback func;	  /* Callback function to be called */
    gpointer data;	  /* User data to pass to the callback function */
  } CallbackData;
  
  
/* Struct to hold a pair of coordinate ranges, one holding ref seq coords, the other match coords */
typedef struct _CoordRange
  {
    int qStart;     /* start coord on the reference (Query) seq */
    int qEnd;       /* end coord on the reference (Query) seq */
    int sStart;     /* start coord on the match (Subject) seq */
    int sEnd;       /* end coord on the match (Subject) seq */
  } CoordRange;
  
  

GdkDrawable*	      widgetGetDrawable(GtkWidget *widget);
void		      widgetSetDrawable(GtkWidget *widget, GdkDrawable *drawable);
gboolean	      widgetGetHidden(GtkWidget *widget);
void		      widgetSetHidden(GtkWidget *widget, const gboolean hidden);
void		      hideUserHiddenWidget(GtkWidget *widget, gpointer data);
void		      widgetClearCachedDrawable(GtkWidget *widget);

gboolean	      onExposePrintableLabel(GtkWidget *label, GdkEventExpose *event, gpointer data);
GtkWidget*	      createLabel(const char *text, const gdouble xalign, const gdouble yalign, const gboolean enableCopyPaste, const gboolean showWhenPrinting);
GdkDrawable*	      createBlankPixmap(GtkWidget *widget);

void		      sortValues(int *val1, int *val2, gboolean forwards);
int		      numDigitsInInt(int val);
gboolean	      parseBlxColor(const char *color, GdkColor *result);
void		      getSelectionColor(GdkColor *origColor, GdkColor *result);
void		      getDropShadowColor(GdkColor *origColor, GdkColor *result);
void		      adjustColorBrightness(GdkColor *origColor, const double factor, GdkColor *result);

gboolean	      mspIsExon(const MSP const *msp);
gboolean	      mspIsIntron(const MSP const *msp);
gboolean	      mspIsSnp(const MSP const *msp);
gboolean	      mspIsBlastMatch(const MSP const *msp);
gboolean	      mspIsCds(const MSP const *msp, const BlxSequence *blxSeq);

gboolean              mspHasSName(const MSP const *msp);
gboolean              mspHasSSeq(const MSP  const *msp);
gboolean              mspHasSCoords(const MSP const *msp);
gboolean              mspHasSStrand(const MSP const *msp);

void		      getCoordRangeExtents(CoordRange *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax);

int		      getRangeLength(const IntRange const *range);
int		      getRangeCentre(const IntRange const *range);
gboolean	      valueWithinRange(const int value, const IntRange const *range);
void		      boundsLimitValue(int *value, const IntRange const *range);
char		      convertBaseToCorrectCase(const char charToConvert, const BlxSeqType seqType);

int		      convertDisplayIdxToDnaIdx(const int inputIdx, 
						const BlxSeqType inputIdxType,
						const int frame, 
						const int baseNum, 
						const int numFrames,
						const gboolean displayRev,
						const IntRange const *dnaIdxRange);

int		      convertDnaIdxToDisplayIdx(const int dnaIdx, 
						const BlxSeqType displaySeqType,
						const int frame,
						const int numFrames, 
						const gboolean displayRev,
						const IntRange const *dnaIdxRange,
						int *baseNum);

int		      mspGetRefFrame(const MSP const *msp, const BlxSeqType seqType);
BlxStrand	      mspGetRefStrand(const MSP const *msp);
BlxStrand	      mspGetMatchStrand(const MSP const *msp);
const char*           mspGetMatchSeq(const MSP const *msp);
char*		      mspGetSeqName(const MSP *msp);
const IntRange const* mspGetRefCoords(const MSP const *msp);
const IntRange const* mspGetMatchCoords(const MSP const *msp);
int		      mspGetQStart(const MSP const *msp);
int		      mspGetQEnd(const MSP const *msp);
int		      mspGetSStart(const MSP const *msp);
int		      mspGetSEnd(const MSP const *msp);
int		      mspGetQRangeLen(const MSP const *msp);
int		      mspGetSRangeLen(const MSP const *msp);
int		      mspGetMatchSeqLen(const MSP const *msp);

char                  getStrandAsChar(const BlxStrand strand);

BlxSequence*          createEmptyBlxSequence(char *fullName);
BlxSequence*          findBlxSequence(GList *seqList, const char *reqdName, const BlxStrand reqdStrand);

int                   roundNearest(const double val);
int		      roundToValue(const int inputVal, const int roundTo);

char		      getRefSeqBase(char *refSeq, 
				    const int qIdx, 
				    const gboolean complement, 
				    const IntRange const *refSeqRange,
				    const BlxSeqType seqType);

int		      getStartDnaCoord(const IntRange const *displayRange, 
				       const int frame,
				       const BlxSeqType displaySeqType, 
				       const gboolean displayRev, 
				       const int numFrames,
				       const IntRange const *refSeqRange);

int		      getEndDnaCoord(const IntRange const *displayRange, 
				     const int frame,
				     const BlxSeqType displaySeqType, 
				     const gboolean displayRev, 
				     const int numFrames,
				     const IntRange const *refSeqRange);

int		      gapCoord(const MSP *msp, 
			       const int qIdx, 
			       const int numFrames, 
			       const BlxStrand strand, 
			       const gboolean displayRev,
			       const gboolean showUnalignedSeq,
			       const gboolean limitUnalignedBases,
			       const gboolean numUnalignedBases);

int		      wildcardSearch(const char *textToSearch, const char *searchStr);

char*		      convertIntToString(const int value);
int		      convertStringToInt(const char *inputStr);
char*		      abbreviateText(const char *inputStr, const int maxLen);
gboolean              stringsEqual(const char *str1, const char *str2, const gboolean caseSensitive);

void		      showMessageDialog(const char *title,  
					const char *messageText,
					GtkWidget *parent,
					const int initWidth,
					const int maxHeight,
					const gboolean wrapText,
					PangoFontDescription *fontDesc);

GtkWidget*		createScrollableTextView(const char *messageText,
						 const gboolean wrapText,
						 PangoFontDescription *fontDesc,
						 int *height);
				    
void		      widgetSetCallbackData(GtkWidget *widget, GtkCallback callbackFunc, gpointer callbackData);
void		      widgetCallCallback(GtkWidget *widget);
void		      widgetCallAllCallbacks(GtkWidget *widget, gpointer data);
void		      onResponseDialog(GtkDialog *dialog, gint responseId, gpointer data);

void		      setPrimaryClipboardText(const char *text);
void		      setDefaultClipboardText(const char *text);
void		      requestPrimaryClipboardText(GtkClipboardTextReceivedFunc callback, gpointer data);
void		      requestDefaultClipboardText(GtkClipboardTextReceivedFunc callback, gpointer data);

int		      parseMatchLine(const char *inputText,
				     char **matchNameOut,
				     int *matchStartOut, 
				     int *matchEndOut, 
				     int *matchLenOut);

GList*		      parseMatchList(const char *inputText);

gint		      runConfirmationBox(GtkWidget *blxWindow, char *title, char *messageText);

const char*	      getSeqVariantName(const char *longName);
char*		      getSeqShortName(const char *longName);
const char*	      sequenceGetFullName(const BlxSequence *seq);
const char*	      sequenceGetVariantName(const BlxSequence *seq);
const char*	      sequenceGetDisplayName(const BlxSequence *seq);
const char*	      sequenceGetShortName(const BlxSequence *seq);
const IntRange const* sequenceGetRange(const BlxSequence *seq);
int		      sequenceGetLength(const BlxSequence *seq);

void		      destroyBlxSequence(BlxSequence *seq);

BlxColor*	      getBlxColor(GList *colorList, const BlxColorId colorId);

void		      defaultMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);
void		      popupMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);

void		      prefixError(GError *error, char *prefixStr, ...);
void                  reportAndClearIfError(GError **error, GLogLevelFlags log_level);

void		      intrangeSetValues(IntRange *range, const int val1, const int val2);

#ifdef DEBUG
void		      debugLogLevel(const int increaseAmt);
#endif

#endif /* _utilities_h_included_ */
