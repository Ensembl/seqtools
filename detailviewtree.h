/*
 *  detailviewtree.h
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#ifndef _detail_view_tree_included_
#define _detail_view_tree_included_

#include <gtk/gtk.h>
#include <SeqTools/blixem_.h>
#include <SeqTools/sequencecellrenderer.h>

#define DETAIL_VIEW_TREE_NAME		  "DetailViewTreeName"
#define DETAIL_VIEW_TREE_CONTAINER_NAME	  "DetailViewTreeContainerName"


/* This struct holds info about a tree header widget. */
typedef struct _TreeColumnHeaderInfo
  {
    GtkWidget *headerWidget;	/* The actual header widget */
    GtkWidget *tree;		/* The tree view that this header belongs to */
    GList *columnIds;		/* A list of columns spanned by this header. (Columns must be adjacent to each other.) */
    GtkCallback refreshFunc;    /* Function to be called when a refresh is requested */
  } TreeColumnHeaderInfo;



typedef struct _TreeProperties 
  {
    GtkWidget *grid;		    /* The grid that this tree corresponds to */
    GtkWidget *detailView;	    /* The detail view that this tree belongs to */

    int readingFrame;		    /* Which reading frame this tree displays */
    GList *treeColumnHeaderList;    /* List of info about the tree column headers */
    
    GHashTable *seqTable;	    /* Hash table to group this tree's MSPs by sequence name. */
    
    GtkTreeModel *mspTreeModel;	    /* Default tree data store, in which each MSP has its own row */
    GtkTreeModel *seqTreeModel;     /* Condensed tree data store, in which multiple MSPs on the same sequence appear in the same row */
  } TreeProperties;
GdkColor exonColour;


/* Public function declarations */
TreeProperties*	  treeGetProperties(GtkWidget *widget);
GtkAdjustment*	  treeGetAdjustment(GtkWidget *tree);
GtkCellRenderer*  treeGetRenderer(GtkWidget *tree);
int		  treeGetCharWidth(GtkWidget *tree);
GtkWidget*	  treeGetGrid(GtkWidget *tree);
Strand		  treeGetStrand(GtkWidget *tree);
gboolean	  treeGetStrandsToggled(GtkWidget *tree);
int		  treeGetNumReadingFrames(GtkWidget *tree);
int		  treeGetSelectedBaseIdx(GtkWidget *tree);
char*		  treeGetRefSeq(GtkWidget *tree);
IntRange*	  treeGetDisplayRange(GtkWidget *tree);
GList*		  treeGetMsps(GtkTreeModel *model, GtkTreeIter *iter);
GtkTreeModel*	  treeGetVisibleDataModel(GtkTreeView *tree);
GtkTreeModel*	  treeGetBaseDataModel(GtkTreeView *tree);
GtkWidget*	  treeGetMainWindow(GtkWidget *tree);
GtkWidget*	  treeGetDetailView(GtkWidget *tree);
int		  treeGetCellXPadding(GtkWidget *tree);
int		  treeGetCellYPadding(GtkWidget *tree);
int		  treeGetCharWidth(GtkWidget *tree);
int		  treeGetCharHeight(GtkWidget *tree);
BlxSeqType	  treeGetSeqType(GtkWidget *tree);
BlxBlastMode	  treeGetBlastMode(GtkWidget *tree);
GHashTable*	  treeGetSeqTable(GtkWidget *tree);

GdkColor*	  treeGetRefSeqColour(GtkWidget *tree, const gboolean selected);
GdkColor*	  treeGetMatchColour(GtkWidget *tree, const gboolean selected);
GdkColor*	  treeGetConsColour(GtkWidget *tree, const gboolean selected);
GdkColor*	  treeGetMismatchColour(GtkWidget *tree, const gboolean selected);
GdkColor*	  treeGetExonColour(GtkWidget *tree, const gboolean selected);
GdkColor*	  treeGetInsertionColour(GtkWidget *tree, const gboolean selected);
GdkColor*	  treeGetExonBoundaryColour(GtkWidget *tree, const gboolean isStart);
int		  treeGetExonBoundaryWidth(GtkWidget *tree);
GdkLineStyle	  treeGetExonBoundaryStyle(GtkWidget *tree, const gboolean isStart);

void		  callFuncOnAllDetailViewTrees(GtkWidget *widget, gpointer data);

void		  deselectAllSiblingTrees(GtkWidget *tree, gboolean includeCurrent);
void		  deselectAllRows(GtkWidget *tree, gpointer data);

void		  treeSortByName(GtkWidget *tree, gpointer data);
void		  treeSortById(GtkWidget *tree, gpointer data);
void		  treeSortByScore(GtkWidget *tree, gpointer data);
void		  treeSortByPos(GtkWidget *tree, gpointer data);
void		  treeSortByGroupOrder(GtkWidget *tree, gpointer data);
void		  refilterTree(GtkWidget *tree, gpointer data);
void		  resortTree(GtkWidget *tree, gpointer data);
void		  refreshTreeHeaders(GtkWidget *tree, gpointer data);
void		  resizeTreeColumns(GtkWidget *tree, gpointer data);
void		  treeUpdateFontSize(GtkWidget *tree, gpointer data);

void		  treeSquashMatches(GtkWidget *tree, gpointer data);
void		  treeUnsquashMatches(GtkWidget *tree, gpointer data);
gboolean	  treeGetMatchesSquashed(GtkWidget *tree);

void		  treeScrollSelectionIntoView(GtkWidget *tree, gpointer data);
void		  selectRowsForSelectedSeqs(GtkWidget *tree, gpointer data);
gboolean	  treeIsSeqSelected(GtkWidget *tree, const char *seqName);

void		  addMspToTree(GtkWidget *tree, MSP *msp);
void		  addSequencesToTree(GtkWidget *tree, gpointer data);

GtkWidget*	  createDetailViewTree(GtkWidget *grid, 
				       GtkWidget *detailView,
				       GtkCellRenderer *renderer,
				       GList **treeList,
				       GList *columnList,
				       BlxSeqType seqType,
				       const char const *refSeqName,
				       const int frame);

void		   treeCreateBaseDataModel(GtkWidget *tree, gpointer data);
void		   treeCreateFilteredDataModel(GtkWidget *tree, gpointer data);

#endif /* _detail_view_tree_included_ */
