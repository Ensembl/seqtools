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

#define UNSET_INT                     -1   /* this value indicates an unset integer */
#define DEFAULT_LABEL_X_PAD           0    /* default x padding to use for header labels */


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



/* Color strings that can be passed to create a GdkColor */
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
#define BLX_VERY_DARK_GREY    "#5f5f5f"
#define BLX_LIGHT_GREY	      "#d3d3d3"
#define BLX_BLACK	      "#000000"
#define BLX_WHITE	      "#ffffff"
#define BLX_RED		      "#ff0000"
#define BLX_LIGHTER_RED	      "#ff5353"
#define BLX_LIGHT_RED	      "#ff7373"
#define BLX_DARK_RED	      "#800000"
#define BLX_VERY_DARK_RED     "#400000"
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


/* This enum contains a list of all the boolean options that the user can toggle on/off */
typedef enum
  {
    BLXFLAG_MIN,		    /* Start index for looping through flags */
  
    BLXFLAG_SQUASH_MATCHES,	    /* Puts all MSPs from the same sequence on the same row in the detail view */
    BLXFLAG_INVERT_SORT,	    /* Inverts the default sort order */
    BLXFLAG_HIGHLIGHT_DIFFS,	    /* Hides matching bases and highlights mis-matching ones */
    BLXFLAG_SHOW_SNP_TRACK,	    /* Shows the SNP track */
    BLXFLAG_SHOW_UNALIGNED,	    /* Shows additional bits of the match sequence that are not part of the aligned section */
    BLXFLAG_SHOW_UNALIGNED_SELECTED,/* Only show unaligned bits of sequence for the currently-selected sequence(s) */
    BLXFLAG_LIMIT_UNALIGNED_BASES,  /* If the show-unaligned-sequence option is on, limits how many bases from the unaligned sequence are shown */
    BLXFLAG_SHOW_POLYA_SITE,        /* Show polyA tails */
    BLXFLAG_SHOW_POLYA_SITE_SELECTED,/* Only show polyA tails for the currently-selected sequence(s) */
    BLXFLAG_SHOW_POLYA_SIG,         /* Show polyA signals in the reference sequence */
    BLXFLAG_SHOW_POLYA_SIG_SELECTED,/* Only show polyA signals for the currently-selected sequence(s) */
    BLXFLAG_SHOW_SPLICE_SITES,	    /* Highlights splice sites in the reference sequence for the currently-selected MSPs */
    BLXFLAG_EMBL_DATA_LOADED,       /* Gets set to true if the full EMBL data is parsed and populated in the MSPs */
    BLXFLAG_SHOW_CDS,               /* True if CDS/UTR regions should be shown; false if plain exons should be shown */
    
    BLXFLAG_NUM_FLAGS		    /* Number of flags, for looping through flags or creating an array */
  } BlxFlag;



/* Function pointer for callback functions used by widgets on dialog boxes. */
typedef gboolean (*BlxResponseCallback)(GtkWidget *widget, const gint responseId, gpointer data);

/* This struct holds generic callback information. It can be stored as a widget
 * property and called on the widget on request (e.g. by a dialog when it applies changes). */
typedef struct _CallbackData
  {
    BlxResponseCallback func;	  /* Callback function to be called */
    gpointer data;                /* User data to pass to the callback function */
  } CallbackData;

/* Custom dialog response types, which can be used in addition to the default types specified by GtkResponseType */
typedef enum
  {
    BLX_RESPONSE_FORWARD, 
    BLX_RESPONSE_BACK
  } BlxResponseType;


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
void		      widgetClearCachedDrawable(GtkWidget *widget, gpointer data);

gboolean	      onExposePrintableLabel(GtkWidget *label, GdkEventExpose *event, gpointer data);
GtkWidget*	      createLabel(const char *text, const gdouble xalign, const gdouble yalign, const gboolean enableCopyPaste, const gboolean showWhenPrinting);
GdkDrawable*	      createBlankPixmap(GtkWidget *widget);

void		      sortValues(int *val1, int *val2, gboolean forwards);
int		      numDigitsInInt(int val);
gboolean              getColorFromString(const char *colorStr, GdkColor *color, GError **error);
void		      getSelectionColor(GdkColor *origColor, GdkColor *result);
void		      getDropShadowColor(GdkColor *origColor, GdkColor *result);
void		      convertToGrayscale(GdkColor *origColor, GdkColor *result);
void		      adjustColorBrightness(GdkColor *origColor, const double factor, GdkColor *result);

gboolean              mspLayerIsVisible(const MSP const *msp);
gboolean              typeIsExon(const BlxMspType mspType);
gboolean              typeIsIntron(const BlxMspType mspType);
gboolean              typeIsMatch(const BlxMspType mspType);
gboolean              typeIsVariation(const BlxMspType mspType);
gboolean	      mspIsExon(const MSP const *msp);
gboolean	      mspIsIntron(const MSP const *msp);
gboolean	      mspIsVariation(const MSP const *msp);
gboolean	      mspIsZeroLenVariation(const MSP const *msp);
gboolean	      mspIsPolyASite(const MSP const *msp);
gboolean	      mspIsPolyASignal(const MSP const *msp);
gboolean	      mspIsBlastMatch(const MSP const *msp);

gboolean              mspHasSName(const MSP const *msp);
gboolean              mspHasSSeq(const MSP  const *msp);
gboolean              mspHasSCoords(const MSP const *msp);
gboolean              mspHasSStrand(const MSP const *msp);
gboolean              mspHasPolyATail(const MSP const *msp, const MSP const *mspList);
gboolean              mspCoordInPolyATail(const int coord, const MSP const *msp, const MSP const *mspList);

void		      getCoordRangeExtents(CoordRange *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax);

int		      getRangeLength(const IntRange const *range);
int		      getRangeCentre(const IntRange const *range);
gboolean	      valueWithinRange(const int value, const IntRange const *range);
gboolean              rangesOverlap(const IntRange const *range1, const IntRange const *range2);
void		      boundsLimitValue(int *value, const IntRange const *range);
void                  boundsLimitRange(IntRange *range, const IntRange const *limit, const gboolean maintainLen);
char		      convertBaseToCorrectCase(const char charToConvert, const BlxSeqType seqType);

void                  convertDisplayRangeToDnaRange(const IntRange const * displayRange, 
                                                    const BlxSeqType displaySeqType,
                                                    const int numFrames,
                                                    const gboolean displayRev,
                                                    const IntRange const *refSeqRange,
                                                    IntRange *result);

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
const char*	      mspGetSName(const MSP *msp);
char*		      blxSequenceGetSummaryInfo(const BlxSequence const *blxSeq);
char*                 mspGetExonTranscriptName(const MSP *msp);
const IntRange const* mspGetRefCoords(const MSP const *msp);
const IntRange const* mspGetMatchCoords(const MSP const *msp);
int		      mspGetQStart(const MSP const *msp);
int		      mspGetQEnd(const MSP const *msp);
int		      mspGetSStart(const MSP const *msp);
int		      mspGetSEnd(const MSP const *msp);
int		      mspGetQRangeLen(const MSP const *msp);
int		      mspGetSRangeLen(const MSP const *msp);
int		      mspGetMatchSeqLen(const MSP const *msp);
const GdkColor*       mspGetColor(const MSP const *msp, GArray *defaultColors, const BlxSequence *blxSeq, const gboolean selected, const gboolean usePrintColors, const gboolean fill);
char*                 mspGetOrganism(const MSP const *msp);
char*                 mspGetGeneName(const MSP const *msp);
char*                 mspGetTissueType(const MSP const *msp);
char*                 mspGetStrain(const MSP const *msp);
char*                 mspGetCoordsAsString(const MSP const *msp);

void                  mspGetFullSRange(const MSP const *msp, 
                                       const gboolean seqSelected,
                                       const gboolean *flags,
                                       const int numUnalignedBases, 
                                       const MSP const *mspList,
                                       IntRange *sSeqRange);

void                  mspGetFullQRange(const MSP const *msp, 
                                       const gboolean seqSelected,
                                       const gboolean *flags,
                                       const int numUnalignedBases, 
                                       const MSP const *mspList,
                                       const int numFrames, 
                                       IntRange *sSeqRange);

char                  getStrandAsChar(const BlxStrand strand);

BlxSequence*          createEmptyBlxSequence(const char *fullName, const char *idTag, GError **error);
BlxSequence*          findBlxSequence(GList *seqList, const char *reqdName, const char *reqdIdTag, const BlxStrand reqdStrand);

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
                               const gboolean seqSelected,
                               const int numUnalignedBases,
                               gboolean *flags,
                               const MSP const *mspList);

int		      wildcardSearch(const char *textToSearch, const char *searchStr);

char*		      convertIntToString(const int value);
char*                 convertDoubleToString(const gdouble value);
int		      convertStringToInt(const char *inputStr);
gboolean	      isWhitespaceChar(const char curChar);
char*		      abbreviateText(const char *inputStr, const int maxLen);
gboolean              stringsEqual(const char *str1, const char *str2, const gboolean caseSensitive);
gboolean	      isValidIupacChar(const char inputChar, const BlxSeqType seqType);

GtkWidget*	      showMessageDialog(const char *title,  
					const char *messageText,
					GtkWidget *parent,
					const int initWidth,
					const int maxHeight,
					const gboolean wrapText,
                                        const gboolean useMarkup,
					PangoFontDescription *fontDesc,
                                        GtkTextView **textView);

void			destroyMessageList();

GtkWidget*		createScrollableTextView(const char *messageText,
						 const gboolean wrapText,
						 PangoFontDescription *fontDesc,
                                                 const gboolean useMarkup,
						 int *height,
                                                 GtkTextView **textViewOut);
				    
void		      widgetSetCallbackData(GtkWidget *widget, BlxResponseCallback callbackFunc, gpointer callbackData);
gboolean              widgetCallAllCallbacks(GtkWidget *widget, gpointer data);
void		      onResponseDialog(GtkDialog *dialog, gint responseId, gpointer data);
void                  onCloseDialog(GtkDialog *dialog, gpointer data);
void                  dialogClearContentArea(GtkDialog *dialog);

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

GtkWidget*	      getNamedChildWidget(GtkWidget *widget, const gchar *searchName);

gint		      runConfirmationBox(GtkWidget *blxWindow, char *title, char *messageText);

const char*	      getSeqVariantName(const char *longName);
char*		      getSeqShortName(const char *longName);
void		      blxSequenceSetName(BlxSequence *seq, const char *fullName);
const char*	      blxSequenceGetFullName(const BlxSequence *seq);
const char*	      blxSequenceGetVariantName(const BlxSequence *seq);
const char*	      blxSequenceGetDisplayName(const BlxSequence *seq);
const char*	      blxSequenceGetShortName(const BlxSequence *seq);
int		      blxSequenceGetLength(const BlxSequence *seq);
char*                 blxSequenceGetSeq(const BlxSequence *seq);
BlxSequence*          blxSequenceGetVariantParent(const BlxSequence *variant, GList *allSeqs);
char*                 blxSequenceGetInfo(BlxSequence *blxSeq, const gboolean allowNewlines, const gboolean dataLoaded);
int		      blxSequenceGetStart(const BlxSequence *seq);
int		      blxSequenceGetEnd(const BlxSequence *seq);

void		      destroyBlxSequence(BlxSequence *seq);

void		      prefixError(GError *error, char *prefixStr, ...);
void                  postfixError(GError *error, char *formatStr, ...);
void                  reportAndClearIfError(GError **error, GLogLevelFlags log_level);

void		      intrangeSetValues(IntRange *range, const int val1, const int val2);

#ifdef DEBUG
void		      debugLogLevel(const int increaseAmt);
#endif

void                  drawHighlightBox(GdkDrawable *drawable, const GdkRectangle const *rect, const gint minWidth, GdkColor *color);

char*                 blxprintf(char *formatStr, ...);
void                  setStatusBarShadowStyle(GtkWidget *statusBar, const char *shadowStyle);


void    gtk_text_buffer_insert_markup             (GtkTextBuffer *buffer,
                                                   GtkTextIter   *iter,
                                                   const gchar   *markup);

void    gtk_text_buffer_insert_markup_with_tag    (GtkTextBuffer *buffer,
                                                   GtkTextIter   *iter,
                                                   const gchar   *markup,
                                                   GtkTextTag    *tag);

void    gtk_text_buffer_set_markup_with_tag       (GtkTextBuffer *buffer,
                                                   const gchar   *markup,
                                                   GtkTextTag    *tag);

void    gtk_text_buffer_set_markup                (GtkTextBuffer *buffer,
                                                   const gchar   *markup);

#endif /* _utilities_h_included_ */
