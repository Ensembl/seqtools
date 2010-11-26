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
#include <seqtoolsUtils/utilities.h>
#include <blixemApp/sequencecellrenderer.h>
#include <blixemApp/detailview.h>

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
    GtkWidget *treeHeader;	    /* The container that contains all the widgets for the tree header */
    GList *treeColumnHeaderList;    /* List of info about the tree column headers */
    gboolean hasSnpHeader;	    /* Whether a SNP track is shown above this tree */
    
    GtkTreeModel *mspTreeModel;	    /* Default tree data store, in which each MSP has its own row */
    GtkTreeModel *seqTreeModel;     /* Condensed tree data store, in which multiple MSPs on the same sequence appear in the same row */
  } TreeProperties;


/* Public function declarations */
BlxViewContext*	  treeGetContext(GtkWidget *tree);
TreeProperties*	  treeGetProperties(GtkWidget *widget);
BlxStrand	  treeGetStrand(GtkWidget *tree);
GList*		  treeGetMsps(GtkTreeModel *model, GtkTreeIter *iter);
GtkTreeModel*	  treeGetBaseDataModel(GtkTreeView *tree);
GtkWidget*	  treeGetBlxWindow(GtkWidget *tree);
int		  treeGetCellXPadding(GtkWidget *tree);
int		  treeGetCellYPadding(GtkWidget *tree);
int               treeGetFrame(GtkWidget *tree);

void		  callFuncOnAllDetailViewTrees(GtkWidget *widget, GtkCallback func, gpointer data);

void              treeSetSortColumn(GtkWidget *tree, gpointer data);
void		  refilterTree(GtkWidget *tree, gpointer data);
void		  resortTree(GtkWidget *tree, gpointer data);
void		  refreshTreeHeaders(GtkWidget *tree, gpointer data);
void		  resizeTreeColumns(GtkWidget *tree, gpointer data);
void		  treeUpdateFontSize(GtkWidget *tree, gpointer data);

void		  treeSquashMatches(GtkWidget *tree, gpointer data);
void		  treeUnsquashMatches(GtkWidget *tree, gpointer data);
gboolean	  treeGetMatchesSquashed(GtkWidget *tree);

gboolean	  treeMoveRowSelection(GtkWidget *tree, const gboolean moveUp, const gboolean shiftModifier);
void		  treeScrollSelectionIntoView(GtkWidget *tree, gpointer data);

void		  addMspToTree(GtkWidget *tree, MSP *msp);
void		  addSequencesToTree(GtkWidget *tree, gpointer data);

GtkWidget*	  createDetailViewTree(GtkWidget *grid, 
				       GtkWidget *detailView,
				       GtkCellRenderer *renderer,
				       GList **treeList,
				       GList *columnList,
				       BlxSeqType seqType,
				       const char const *refSeqName,
				       const int frame,
				       const gboolean includeSnpTrack);

void		   treeCreateBaseDataModel(GtkWidget *tree, gpointer data);
void		   treeCreateFilteredDataModel(GtkWidget *tree, gpointer data);

#endif /* _detail_view_tree_included_ */
