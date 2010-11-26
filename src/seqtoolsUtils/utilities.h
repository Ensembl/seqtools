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

#define UNSET_INT                     -1   /* this value indicates an unset integer */
#define DEFAULT_LABEL_X_PAD           0    /* default x padding to use for header labels */

/* Really the buffers that use this should be dynamic but I'm not going to do that, this
 * code is so poor that it doesn't warrant the effort.... */
#define NAMESIZE      12
#define MAXLINE       10000


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

/* Generic SeqTools error domain */
#define SEQTOOLS_ERROR g_quark_from_string("SeqTools")

/* Error codes */
typedef enum
{
  SEQTOOLS_ERROR_PARSING_COLOR,	      /* error parsing color string */
  SEQTOOLS_ERROR_SEQ_SEGMENT,	      /* error getting the requested segment of a sequence */
  SEQTOOLS_ERROR_NO_STYLE	      /* style does not exist */
} SeqToolsError;


#if defined(MACINTOSH)
#define      SYSERR_FORMAT             "system error %d"
#else
#define      SYSERR_FORMAT             "system error %d - %s"
#endif


/* Special characters for displaying in sequences */
#define SEQUENCE_CHAR_GAP    '.'   /* represents a gap in the match sequence */
#define SEQUENCE_CHAR_PAD    '-'   /* used for padding when the sequence is unavailable */
#define SEQUENCE_CHAR_BLANK  '-'   /* used to display a blank when we're not interested in what the actual base is */
#define SEQUENCE_CHAR_STOP   '*'   /* STOP codon */
#define SEQUENCE_CHAR_MET    'M'   /* MET codon */


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


typedef struct _BlxColor
  {
    char *name;			  /* meaningful name for the color e.g. "Match" */
    char *desc;			  /* meaningful description for what the color is used for e.g. "background color for exact matches" */
    gboolean transparent;	  /* if this is true, the colors below are not specified and the background color should be used instead */
    
    GdkColor normal;		  /* the color in normal operation */
    GdkColor selected;		  /* the color in a selected state */
    GdkColor print;		  /* the color used for printing */
    GdkColor printSelected;	  /* the selected-state color used for printing */
  } BlxColor;
  
  
/* This handle holds a list of pointers to all memory allocated via this handle. Use handleDestroy
 * to free the handle and all its allocated memory. */
typedef GSList* BlxHandle;
  
  
/* This struct is used to pass user data to the message handlers */
typedef struct _BlxMessageData
{
  GtkWindow *parent;
  GtkStatusbar *statusBar;
} BlxMessageData;
  
  
/* Define a drawing style for an MSP */
typedef struct _BlxStyle
  {
    char *styleName;
    BlxColor fillColor;
    BlxColor lineColor;
    BlxColor fillColorUtr; /* used only for transcript features; fillColor and lineColor are for CDS portions and these are for UTR portions */
    BlxColor lineColorUtr; /* used only for transcript features; fillColor and lineColor are for CDS portions and these are for UTR portions */
  } BlxStyle;

/* Utility structs to hold a range of integers/doubles */
typedef struct _IntRange
  {
    int min;
    int max;
  } IntRange ;

typedef struct _DoubleRange
  {
    gdouble min;
    gdouble max;
  } DoubleRange ;
  

/* Fundamental strand direction. */
typedef enum
  {
    BLXSTRAND_NONE, 
    BLXSTRAND_FORWARD, 
    BLXSTRAND_REVERSE
  } BlxStrand ;
  
/* Fundamental type of sequence (DNA really means nucleotide, because it could be RNA as well). */
typedef enum
  {
    BLXSEQ_INVALID, 
    BLXSEQ_DNA, 
    BLXSEQ_PEPTIDE
  } BlxSeqType ;

/* Fundamental blast mode used */
typedef enum
  {
    BLXMODE_UNSET, 
    BLXMODE_BLASTX, 
    BLXMODE_TBLASTX, 
    BLXMODE_BLASTN, 
    BLXMODE_TBLASTN, 
    BLXMODE_BLASTP
  } BlxBlastMode ;
  

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

BlxSeqType            determineSeqType(char *seq);
void                  argvAdd(int *argc, char ***argv, char *s);
char*                 getSystemErrorText();
gpointer              handleAlloc(BlxHandle *handle, size_t numBytes);
BlxHandle             handleCreate();
void                  handleDestroy(BlxHandle *handle);
BlxStyle*             getBlxStyle(const char *styleName, GSList *styles, GError **error);

void		      sortValues(int *val1, int *val2, gboolean forwards);
int		      numDigitsInInt(int val);
gboolean              getColorFromString(const char *colorStr, GdkColor *color, GError **error);
void		      getSelectionColor(GdkColor *origColor, GdkColor *result);
void		      getDropShadowColor(GdkColor *origColor, GdkColor *result);
void		      convertToGrayscale(GdkColor *origColor, GdkColor *result);
void		      adjustColorBrightness(GdkColor *origColor, const double factor, GdkColor *result);

void		      getCoordRangeExtents(CoordRange *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax);

int		      getRangeLength(const IntRange const *range);
int		      getRangeCentre(const IntRange const *range);
void                  centreRangeOnCoord(IntRange *range, const int coord, const int length);
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

char                  getStrandAsChar(const BlxStrand strand);

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

int		      wildcardSearch(const char *textToSearch, const char *searchStr);

char*		      convertIntToString(const int value);
char*                 convertDoubleToString(const gdouble value, const int numDp);
int		      convertStringToInt(const char *inputStr);
gboolean	      isWhitespaceChar(const char curChar);
char*		      abbreviateText(const char *inputStr, const int maxLen);
gboolean              stringsEqual(const char *str1, const char *str2, const gboolean caseSensitive);
gboolean	      isValidIupacChar(const char inputChar, const BlxSeqType seqType);
void                  stringProtect(FILE *file, const char *string);
char*                 stringUnprotect(char **textp, char *target);

int                   invertCoord(const int coord, const IntRange const *range, const gboolean invert);

void                  popupMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);
void		      defaultMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);

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
GtkWidget*            getPersistentDialog(GtkWidget* dialogList[], const int dialogId);
void                  addPersistentDialog(GtkWidget* dialogList[], const int dialogId, GtkWidget *widget);

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

BlxColor*	      getBlxColor(GArray *defaultColors, const int colorId);
GdkColor*	      getGdkColor(int colorId, GArray *defaultColors, const gboolean selected, const gboolean usePrintColors);
const GdkColor*	      blxColorGetColor(const BlxColor *blxColor, const gboolean selected, const gboolean usePrintColors);
char*		      convertColorToString(GdkColor *color);
void		      destroyBlxColor(BlxColor *blxColor);

void		      createBlxColor(GArray *defaultColors,
				     int colorId,
				     const char *name, 
				     const char *desc, 
				     const char *normalCol, 
				     const char *printCol,
				     const char *normalColSelected,
				     const char *printColSelected);

gdouble		      pixelsPerBase(const gint displayWidth, 
				    const IntRange const *displayRange);

gint		      convertBaseIdxToRectPos(const gint dnaIdx, 
					      const GdkRectangle const *rect, 
					      const IntRange const *dnaDispRange,
					      const gboolean displayRev,
					      const gboolean clip);

gchar*			  getSequenceSegment(const char const *dnaSequence,
                                             IntRange *qRange,
					     const BlxStrand strand,
					     const BlxSeqType srcSeqType,
					     const BlxSeqType destSeqType,
					     const int frame,
					     const int numFrames,
					     const IntRange const *refSeqRange,
					     const BlxBlastMode blastMode,
					     char **geneticCode,
					     const gboolean displayRev,
					     const gboolean reverseResult,
					     const gboolean allowComplement,
					     GError **error);

const char*			   findFixedWidthFont(GtkWidget *widget);
const char*			   findFixedWidthFontFamily(GtkWidget *widget, GList *pref_families);
void                               getFontCharSize(GtkWidget *widget, PangoFontDescription *fontDesc, gdouble *width, gdouble *height);

GtkWidget*                         createEmptyButtonBar(GtkToolbar **toolbar);
void                               makeToolbarButton(GtkToolbar *toolbar, char *label, char *stockId, char *tooltip, GtkSignalFunc callback_func, gpointer data);


/* seqtoolsWebBrowser.c */
gboolean                           seqtoolsLaunchWebBrowser(const char *link, GError **error);


/* translate.c */
char*                              blxTranslate(const char *seq, char **code);
void                               blxComplement(char *seq) ;    
char*                              revComplement(char *comp, char *seq) ;


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
