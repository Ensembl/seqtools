/*  File: belvuTree.c
 *  Author: Gemma Barson, 2011-05-06
 *  Copyright (c) 2011 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * SeqTools is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * ---------------------------------------------------------------------------
 * This file is part of the SeqTools sequence analysis package, 
 * written by
 *      Gemma Barson      (Sanger Institute, UK)  <gb10@sanger.ac.uk>
 * 
 * based on original code by
 *      Erik Sonnhammer   (SBC, Sweden)           <Erik.Sonnhammer@sbc.su.se>
 * 
 * and utilizing code taken from the AceDB and ZMap packages, written by
 *      Richard Durbin    (Sanger Institute, UK)  <rd@sanger.ac.uk>
 *      Jean Thierry-Mieg (CRBM du CNRS, France)  <mieg@kaa.crbm.cnrs-mop.fr>
 *      Ed Griffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: See belvuTree.h
 *----------------------------------------------------------------------------
 */


#include "belvuApp/belvuTree.h"
#include "belvuApp/belvuWindow.h"
#include "seqtoolsUtils/utilities.h"
#include <string.h>
#include <ctype.h> /* for isspace etc. */
#include <math.h>


#define BELVU_TREE_WINDOW_NAME                  "BelvuTreeWindow"
#define DEFAULT_TREE_WINDOW_WIDTH_FRACTION      0.5    /* default width of tree window (as fraction of screen width) */
#define DEFAULT_TREE_WINDOW_HEIGHT_FRACTION     0.4   /* default height of tree window (as fraction of screen height) */
#define DEFAULT_XPAD                            10
#define DEFAULT_YPAD                            10


/*  BLOSUM62 930809
 
 #  Matrix made by matblas from blosum62.iij
 #  * column uses minimum score
 #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
 #  Blocks Database = /data/blocks_5.0/blocks.dat
 #  Cluster Percentage: >= 62
 #  Entropy =   0.6979, Expected =  -0.5209
 
 Note: to use with a2b[], always subtract 1 from the values !!!!
 
 A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X  \* */ 
static int BLOSUM62[24][24] = {
{4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
{0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
{-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
{-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
{ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -1}
};


/* ASCII-to-binary translation table
 Note: to use with BLOSUM62[], always subtract 1 from the values !!!! */
#undef NA
#define NA 0
static int a2b[] =
{
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA, 1,NA, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
15, 6, 2,16,17,NA,20,18,NA,19,NA,NA,NA,NA,NA,NA,
NA, 1,NA, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
15, 6, 2,16,17,NA,20,18,NA,19,NA,NA,NA,NA,NA,NA,

NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA 
};




/* Properties specific to the belvu tree */
typedef struct _BelvuTreeProperties
  {
    BelvuContext *bc;	            /* The belvu context */
    
    GtkWidget *treeArea;            /* Drawing widget for the tree */
    GdkRectangle treeRect;          /* Specifies the actual area in which we'll draw the tree within treeAre */
    
    TreeNode *treeHead;             /* The head node of the tree */
    
    int charWidth;
    int charHeight;
    
    double lineWidth;               /* The line width of the branches */
    double treeScale;               /* The tree scale */
    gboolean showBranchLen;         /* Whether to show the branch lengths on the branches */
    gboolean showOrganism;          /* Whether to show the organism name */
    
  } BelvuTreeProperties;



/* Local function declarations */
static Tree*                        createEmptyTree();

/***********************************************************
 *                         Properties                      *
 ***********************************************************/

static BelvuTreeProperties* belvuTreeGetProperties(GtkWidget *widget)
{
  return widget ? (BelvuTreeProperties*)(g_object_get_data(G_OBJECT(widget), "BelvuTreeProperties")) : NULL;
}

static void onDestroyBelvuTree(GtkWidget *belvuTree)
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
  
  if (properties)
    {
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(belvuTree), "BelvuTreeProperties", NULL);
    }
}


/* Create the properties struct and initialise all values. */
static void belvuTreeCreateProperties(GtkWidget *belvuTree, 
                                      BelvuContext *bc,
                                      GtkWidget *treeArea,
                                      TreeNode *treeHead)
{
  if (belvuTree)
    {
      BelvuTreeProperties *properties = g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->treeArea = treeArea;
      properties->treeHead = treeHead;

      properties->charHeight = 14;
      properties->charWidth = 8;

      properties->lineWidth = bc->treeLineWidth;
      properties->treeScale = bc->treeScale;
      properties->showBranchLen = bc->treeShowBranchlen;
      properties->showOrganism = bc->treeShowOrganism;
      
      g_object_set_data(G_OBJECT(belvuTree), "BelvuTreeProperties", properties);
      g_signal_connect(G_OBJECT(belvuTree), "destroy", G_CALLBACK (onDestroyBelvuTree), NULL);
    }
}


/***********************************************************
 *                       Tree utils                        *
 ***********************************************************/

static void treeTraverse(BelvuContext *bc, TreeNode *node, void (*func)(BelvuContext *bc, TreeNode *treeNode)) 
{
  if (!node) 
    return;
  
  treeTraverse(bc, node->left, func);
  func(bc, node);
  treeTraverse(bc, node->right, func);
}


/***********************************************************
 *                   Tree bootstrapping                    *
 ***********************************************************/

static int BSorder(gconstpointer xIn, gconstpointer yIn)
{
  const BootstrapGroup *x = (const BootstrapGroup*)xIn;
  const BootstrapGroup *y = (const BootstrapGroup*)yIn;
  
  return strcmp(x->s, y->s);
}


static BootstrapGroup* createEmptyBootstrapGroup()
{
  BootstrapGroup *result = g_malloc(sizeof *result);
  
  result->node = NULL;
  result->s = NULL;
  
  return result;
}


/* Combines left and right sequence groups and insert to bootstrap group list */
static GArray* fillBootstrapGroups(BelvuContext *bc, TreeNode *node, TreeNode *rootnode, int maintree) 
{
  GArray *result = NULL;
  
  BootstrapGroup *BS = createEmptyBootstrapGroup();
  
  if (!node->name)
    {
      /* Internal node */
      /* Combine left node sequences into right node array */
      
      GArray *left = fillBootstrapGroups(bc, node->left, rootnode, maintree);
      GArray *right = fillBootstrapGroups(bc, node->right, rootnode, maintree);
      
      /* Nothing to do for root */
      if (node != rootnode)
        {      
          /* Combine left and right groups */
          int i = 0;
          for (i = 0 ; i < left->len; ++i) 
            {
              char *s = g_array_index(left, char*, i);
              g_array_append_val(right, s);
              g_array_sort(right, strcmp_);
            }
          
          g_array_unref(left);
          
          /* Create string with group members */
          int ssize = 0;
          for (i = 0 ; i < right->len ; ++i) 
            ssize += (strlen(g_array_index(right, char*, i)) + 1);
          
          BS->s = g_malloc(ssize+1);
          
          for (i = 0 ; i < right->len ; ++i) 
            {
              strcat(BS->s, g_array_index(right, char*, i));
              strcat(BS->s, " ");
            }
          
          /* printf("   New bootstrap group: %s\n", BS->s); */
          
          if (maintree) 
            {
              /* Associate this node with the group string */
              BS->node =  node;
              
              /* Add group string to array of bootstrap groups */
              g_array_append_val(bc->bootstrapGroups, BS);
              g_array_sort(bc->bootstrapGroups, BSorder);
            }
          else
            {
              /* Find group string and increment counter if exists */
              BootstrapGroup *BS2;
              
              int ip = 0;
              if (arrayFind(bc->bootstrapGroups, BS, &ip, (void *)BSorder)) 
                {
                  BS2 = &g_array_index(bc->bootstrapGroups, BootstrapGroup, ip);
                  BS2->node->boot++;
                  /* printf("Found bootgroup %s\n", BS->s); */
                }
              else
                {
                  /* printf("Did not find bootgroup %s\n", BS->s); */
                }
              
              g_free(BS->s);
              g_free(BS);
            }
          
          result = right;
        }
    }
  else
    {
      /* Leaf node - return as embryonic array with one entry */
      GArray *leaf = g_array_sized_new(FALSE, FALSE, sizeof(char*), bc->alignArr->len);
      g_array_append_val(leaf, node->name);
      g_array_sort(leaf, strcmp_);
      result = leaf;
    }
  
  return result;
}


static void normaliseBootstraps(BelvuContext *bc, TreeNode *node) 
{
  node->boot = 100*node->boot/(double)bc->treebootstraps;
}


/* Bootstrap the internal nodes in a tree.
 1. Set up a list (array) of all internal nodes with leaf content and pointer to node
 2. In bootstrap tree, for each internal node, check if its contents exists in list. If so, increment node's bootstrap count
 3. Turn increments to percentages
 */
static void treeBootstrapStats(BelvuContext *bc, TreeNode *tree)
{
  /* Traverse tree, fill array bootstrapGroups */
  bc->bootstrapGroups = g_array_sized_new(FALSE, FALSE, sizeof(BootstrapGroup), bc->alignArr->len);
  fillBootstrapGroups(bc, tree, tree, 1);
  
  treeBootstrap(bc);
  
  /* Turn increments to percentages */
  treeTraverse(bc, tree, normaliseBootstraps);
}


void treeBootstrap(BelvuContext *bc)
{
  Tree *treeStruct = createEmptyTree();
  
  separateMarkupLines(bc);
  GArray *alignArrTmp = copyAlignArray(bc->alignArr);
  
#if defined(__CYGWIN__) || defined(DARWIN)
  srand(time(0));
#else
  srand48(time(0));
#endif
  
  int iter = 0;
  for (iter = 0; iter < bc->treebootstraps; ++iter) 
    {
      /* Make bootstrapped Alignment */
      int i = 0;
      for (i = 0; i < bc->maxLen; ++i) 
        {
          int src = 0;
          
#if defined(__CYGWIN__) || defined(DARWIN)
          src = (int)((((double)rand())/RAND_MAX) * bc->maxLen);
#else
          src = (int)(drand48() * bc->maxLen);
#endif
          
          columnCopy(bc->alignArr, i, alignArrTmp, src);
        }
      
      
      treeStruct->head = treeMake(bc, 0);
      
      if (bc->outputBootstrapTrees) 
        {
          if (bc->treebootstrapsDisplay) 
            {
              createBelvuTreeWindow(bc, treeStruct->head);
            }
          else
            {
              printf("\n");
              treePrintNH(treeStruct, treeStruct->head, stdout);
              printf(";\n");
            }
        }
      else
        {
          /* Collect bootstrap statistics */
          fillBootstrapGroups(bc, treeStruct->head, treeStruct->head, 0);
        }
    }
  
  /* Restore alignment */
  int i = 0;
  for (i = 0; i < bc->alignArr->len; ++i)
    {
      strcpy(g_array_index(bc->alignArr, ALN, i).seq, g_array_index(alignArrTmp, ALN, i).seq);
    }
}


/***********************************************************
 *                   Business logic                        *
 ***********************************************************/

/*
 * Expect format:
 *
 * name1 name2 name3
 * 1-1   1-2   1-3
 * 2-1   2-2   2-3
 * 3-1   3-2   3-3
 */
static void treeReadDistances(BelvuContext *bc, double **pairmtx)
{
  char   *p;
  int    i = 0, j = 0 ;
  double  d;
  char line[MAXLENGTH+1];
  
  while (!feof (bc->treeReadDistancesPipe))
    { 
      if (!fgets (line, MAXLENGTH, bc->treeReadDistancesPipe))
	break;
      
      j = 0;
      
      while ( (p = strtok(j ? 0 : line, " \t")) && !isspace(*p))
	{
	  d = atof(p);
	  DEBUG_OUT("%d  %d  %f\n", i, j, d);
	  pairmtx[i][j] = d;
	  j++;
	}
      
      if (j != bc->alignArr->len)
	g_error("nseq = %d, but read %d distances in row %d", bc->alignArr->len, j, i);
      
      i++;
    }
  
  if (j != bc->alignArr->len)
    g_error("nseq = %d, but read %d distance rows", bc->alignArr->len, i);
  
  return ;
}


/* Correct an observed distance to an estimated 'true' distance.
 
 Adapted from Lasse Arvestad's lapd program
 
 od = observed distance in percent 
 
 */
static double treeJUKESCANTOR(double od)
{
  double cd;
  
  od /= 100;
  
  od = 20.0/19.0*od;
  if  (od > 0.95) od = 0.95; /* Limit to 300 PAM */
  
  cd = -19.0/20.0 * log(1-od);
  
  return cd*100;
}



/* Correct an observed distance to an estimated 'true' distance.
 
 Adapted from Lasse Arvestad's lapd program
 
 od = observed distance in percent 
 
 */
static double treeKimura(double od)
{
  double adjusted, cd;
  
  od /= 100;
  
  adjusted = od + 0.2 * od*od;
  
  if (adjusted > 0.95) adjusted = 0.95;   /* Limit to 300 PAM */
  
  cd = -log(1 - adjusted);
  
  return cd*100;
}



/* Correct an observed distance to an estimated 'true' distance.
 
 Based on Christian Storm's 5th order polynomial curve fitting of ROSE simulated data
 
 od = observed distance in percent 
 
 */
static double treeSTORMSONN(double od)
{
  double cd;			/* Corrected distance */
  
  if (od > 91.6) return 1000.0;	/* Otherwise log of negative value below */
  
  od /= 100;
  
  cd= -log(1 
           -0.95844*od 
           -0.69957*od*od
           +2.4955*od*od*od
           -4.6353*od*od*od*od
           +2.8076*od*od*od*od*od);
  
  /* printf(" od=%f  cd=%f\n", od, cd); */
  
  if (cd > 3.0) cd=3.0; /* Limit to 300 PAM */
  
  return cd*100;
}


static double treeSCOREDIST(char *seq1, char *seq2, BelvuContext *bc)
{
  int len,
  sc = 0, 
  s1sc = 0, 
  s2sc = 0;
  double od, cd, 
  maxsc, 
  expect; 
  char *s1, *s2;
  
#define mtx BLOSUM62
  
  
  s1 = seq1;
  s2 = seq2;
  
  
  /* Calc scores */
  for (len=0; *s1; s1++, s2++) 
    {
      if (bc->penalize_gaps) 
        {
          if (isGap(*s1) || isGap(*s2)) sc -= 0.6;	
        }
      
      if (!isGap(*s1) && !isGap(*s2)) 
        {
          sc += mtx[a2b[(unsigned char)(*s1)]-1][a2b[(unsigned char)(*s2)]-1];
          s1sc += mtx[a2b[(unsigned char)(*s1)]-1][a2b[(unsigned char)(*s1)]-1];
          s2sc += mtx[a2b[(unsigned char)(*s2)]-1][a2b[(unsigned char)(*s2)]-1];
          len++;
        }
    }
  
  maxsc = (s1sc + s2sc) / 2.0;
  
  /* Calc expected score */
  expect =  -0.5209 * len;
  
  od = ((double)sc - expect) / (maxsc - expect);
  
  if (!len || od < 0.05) od = 0.05;  /* Limit to 300 PAM;  len==0 if no overlap */
  if (od > 1.0) od = 1.0; 
  
  cd = -log(od);
  
  /* printf ("len=%d  sc=%.2f  maxsc=%.2f  expect=%.2f  maxsc-expect=%.2f  od=%.3f\n", len, sc, maxsc, expect, maxsc-expect, od); */
  
  cd = cd * 100.0 * 1.337 /* Magic scaling factor optimized for Dayhoff data */ ;
  
  if (cd > 300) cd = 300;  /* Limit to 300 PAM */
  
  /*    if (!bc->penalize_gaps) {
   g_free(s1);
   g_free(s2);
   } */
  
  return cd;
}



/* Sum branchlengths, allow going up parents
 
 If going up parents, arg1 = parent  ;  arg2 = child.
 
 If you're not interested in going up parents, simply call me with
 the same node in both arguments.
 
 */
static double treeSize3way(TreeNode *node, TreeNode *fromnode) 
{
  int root = 0;
  TreeNode *left = NULL;
  TreeNode *right = NULL;
  double len = 0.0;
  
  if (!node) 
    return 0.0;
  
  /* Get the correct branch length */
  if (node->left == fromnode || node->right == fromnode) /* Going up the tree */
    len = fromnode->branchlen;
  else
    len = node->branchlen;
  
  if (node->left == fromnode) 
    {
      left = node->parent;
      if (!left) 
        root = 1;
    }
  else
    {
      left = node->left;
    }
  
  if (node->right == fromnode) 
    {
      right = node->parent;
      if (!right) 
        root = 1;
    }
  else 
    {
      right = node->right;
    }
  
  if (root)
    {
      double retval;
      
      /* go across tree root */
      if (left) /* Coming from right */
        retval = treeSize3way(left, left);
      else  /* Coming from left */
        retval = treeSize3way(right, right);
      
      DEBUG_OUT("Returning (%.1f + %.1f) = %.1f\n", fromnode->branchlen, retval, fromnode->branchlen + retval);
      
      return fromnode->branchlen + retval;
    }
  else 
    {
      double 
      l = treeSize3way(left, node),
      r = treeSize3way(right, node),
      retval = (l + r)/2.0 + len;
      
      DEBUG_OUT("Returning (%.1f + %.1f)/2 + %.1f = %.1f\n", l, r, len, retval);
      
      return retval;
    }
}


/* Calculate the difference between left and right trees if the tree
 were to be rerooted at this node.
 
 What is calculated (bal) is the difference between 'left' and
 'right' subtrees minus the length of the branch itself. If this
 difference is negative, a perfectly balanced tree can be made.  For
 imperfectly balanced trees we want to root at the branch that gives
 the best balance.  However, perfectly balanced trees are all
 'perfect' so here we chose the branch with most equal subtrees.
 
 Actually it is not "left" and "right" but "down" and "up" subtrees.  */
static void treeCalcBalance(BelvuContext *bc, TreeNode *node) 
{
  double bal, lweight, rweight;
  
  if (node == bc->treeBestBalancedNode) 
    return;
  
  DEBUG_OUT("Left/Downstream weight\n");
  
  lweight = treeSize3way(node, node);
  
  DEBUG_OUT("Right/Upstream weight\n");
  
  rweight = treeSize3way(node->parent, node);
  
  bal = fabsf(lweight - rweight) - node->branchlen;
  
  DEBUG_OUT("Node=%s (branchlen = %.1f).  Weights = %.1f  %.1f. Bal = %.1f\n", 
            node->name, node->branchlen, lweight, rweight, bal);
  
  if (bal < bc->treeBestBalance) 
    { /* better balance */
      if (bc->treeBestBalance > 0.0 ||
          /* If previous tree was not perfectly balanced, or
       If previous tree was perfectly balanced - choose root with best subtree balance */
          fabsf(lweight - rweight) < bc->treeBestBalance_subtrees)
        {
          DEBUG_OUT("            %s has better balance %.1f < %.1f\n", node->name, bal, bc->treeBestBalance);
          
          bc->treeBestBalancedNode = node;
          bc->treeBestBalance = bal;
          bc->treeBestBalance_subtrees = fabsf(lweight - rweight);
        }
    }
}


static TreeNode *treeParent2leaf(TreeNode *newparent, TreeNode *curr)
{
  if (!curr->parent) 
    { /* i.e. the old root */
      if (curr->left == newparent)
        {
          newparent->branchlen += curr->right->branchlen; /* Add second part of this vertex */
          return curr->right;
        }
      else
        {
          newparent->branchlen += curr->left->branchlen;
          return curr->left;
        }
    } 
  else
    {
      if (curr->left == newparent) 
        {
          /* Change the link to the new parent to the old parent */
          curr->left = treeParent2leaf(curr, curr->parent);
          curr->left->branchlen = curr->branchlen;
        }
      else
        {
          curr->right = treeParent2leaf(curr, curr->parent);
          curr->right->branchlen = curr->branchlen;
        }
    }
  
  return curr;
}


/* Balance the two sides of a tree branch by the weights of the two subtrees.
 The branch has two sides: llen and rlen.
 
 Rationale: The correct center point balances the treesizes on both sides.
 
 Method: Calculate half the unbalance, add it to the appropriate side.
 
 No real theory for this, but it seems to work in easy cases 
 */
static void treeBalanceByWeight(TreeNode *lnode, TreeNode *rnode, double *llen, double *rlen)
{
  double adhocRatio = 0.95;
  
  double 
  halfbal, 
  branchlen = *rlen+*llen;
  double lweight = treeSize3way(lnode, lnode) /*- lnode->branchlen*/;
  double rweight = treeSize3way(rnode, rnode) /*- rnode->branchlen*/;
  
  halfbal = fabsf(lweight-rweight) / 2.0;
  
  if (halfbal < *llen && halfbal < *rlen) {
    
    if (lweight < rweight) {
      *llen += halfbal;
      *rlen -= halfbal;
    }
    else {
      *llen -= halfbal;
      *rlen += halfbal;
    }
  }
  else {
    /* The difference is larger than the final branch -
     give nearly all weight to the shorter one (cosmetic hack) */
    if (lweight < rweight) {
      *llen = branchlen*adhocRatio;
      *rlen = branchlen*(1.0 - adhocRatio);
    }
    else {
      *rlen = branchlen*adhocRatio;
      *llen = branchlen*(1.0 - adhocRatio);
    }
  }
}


static void fillParents(TreeNode *parent, TreeNode *node)
{
  if (!node) 
    return;
  
  node->parent = parent;
  fillParents(node, node->left);
  fillParents(node, node->right);
}


static char *fillOrganism(TreeNode *node)
{
  char 
  *leftorganism,
  *rightorganism;
  
  if (node->name) 
    return node->organism;
  
  leftorganism = fillOrganism(node->left);
  rightorganism = fillOrganism(node->right);
  
  /* printf("\nFill: (left=%s, right=%s):  ", leftorganism, rightorganism); */
  
  node->organism = (leftorganism == rightorganism ? leftorganism : NULL);
  
  return node->organism;
}


/* Allocate memory for a new TreeNode and initialise its contents to empty values */
static TreeNode* createEmptyTreeNode()
{
  TreeNode *result = g_malloc(sizeof *result);
  
  result->dist = 0.0;
  result->branchlen = 0.0;
  result->boot = 0.0;
  result->left = NULL;
  result->right = NULL;
  result->parent = NULL;
  result->name = NULL;
  result->organism =NULL;
  result->aln = NULL;
  result->box = 0;
  result->color = 0;
  
  return result;
}


/* Rerooting works roughly this way:
 
 - A new node is created with one child being the node chosen as new root 
 and the other child the previous parent of it.  The sum of left and right
 branch lengths of the new root should equal the branch length of the chosen node.
 - All previous parent nodes are visited and converted so that:
 1. The previous parent node becomes a child.
 2. The new parent is the node that calls.
 3. The branchlength of the previous parent node is assigned to its parent.
 - When the previous root is reached, it is deleted and the other child
 becomes the child of the calling node.
 
 Note that treeReroot destroys the old tree, since it reuses the nodes.  
 Use treecpy if you still need it for something later.
 */
static TreeNode *treeReroot(BelvuContext *bc, TreeNode *node)
{
  TreeNode *newroot = createEmptyTreeNode();
  
  newroot->left = node;
  newroot->right = treeParent2leaf(node, node->parent);
  
  fillParents(newroot, newroot->left);
  fillParents(newroot, newroot->right);
  
  newroot->left->branchlen = newroot->right->branchlen = node->branchlen / 2.0;
  treeBalanceByWeight(newroot->left, newroot->right, 
                      &newroot->left->branchlen, &newroot->right->branchlen);
  
  bc->treeHead = newroot;
  return newroot;
}


/* Find the node which has most equal balance, return new tree with this as root.
 */
static TreeNode *treeFindBalance(BelvuContext *bc, TreeNode *tree) 
{
  double lweight = treeSize3way(tree->left, tree->left);
  double rweight = treeSize3way(tree->right, tree->right);
  
  bc->treeBestBalancedNode = tree;
  bc->treeBestBalance = fabsf(lweight - rweight);
  
  bc->treeBestBalance_subtrees = 
  fabsf((lweight - tree->left->branchlen) - (rweight - tree->right->branchlen));
  
  DEBUG_OUT("Initial weights = %.1f  %.1f. Bal = %.1f\n", lweight, rweight, bc->treeBestBalance);
  
  treeTraverse(bc, tree, treeCalcBalance);
  
  if (bc->treeBestBalancedNode == tree)
    return tree;
  else
    return treeReroot(bc, bc->treeBestBalancedNode);
}


TreeNode *treeMake(BelvuContext *bc, const gboolean doBootstrap)
{
  TreeNode *newnode = NULL ;
  int maxi = -1, maxj = -1;
  double maxid = 0.0, pmaxid, **pairmtx, **Dmtx, **curMtx, *src, *trg, 
  *avgdist,		/* vector r in Durbin et al */
  llen = 0, rlen = 0;
  TreeNode **node ;					    /* Array of (primary) nodes.  Value=0 => stale column */
  static BlxHandle handle = 0;
  BlxHandle treeHandle = NULL ;
  
  if (doBootstrap)
    treeHandle = handle;
  
  /* Allocate memory */
  if (treeHandle)
    handleDestroy(&treeHandle);
  
  treeHandle = handleCreate();
  pairmtx = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double *));
  Dmtx = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double *));
  node = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(TreeNode *));
  avgdist = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double));

  int i = 0;
  for (i = 0; i < bc->alignArr->len; ++i)
    {
      pairmtx[i] = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double));
      Dmtx[i] = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double));
      node[i] = handleAlloc(&treeHandle, sizeof(TreeNode));
      node[i]->name =  handleAlloc(&treeHandle, strlen(g_array_index(bc->alignArr, ALN, i).name)+50);
      
      ALN *aln_i = &g_array_index(bc->alignArr,ALN, i);
      
      if (!bc->treeCoordsOn) 
	{
	  sprintf(node[i]->name, "%s", aln_i->name);
	}
      else
	{
	  sprintf(node[i]->name, "%s/%d-%d", 
		  aln_i->name,
		  aln_i->start,
		  aln_i->end);
	}
      
      node[i]->aln = aln_i;
      node[i]->organism = aln_i->organism;

      node[i]->dist = 0.0;
      node[i]->branchlen = 0.0;
      node[i]->boot = 0.0;
      node[i]->left = NULL;
      node[i]->right = NULL;
      node[i]->parent = NULL;
      node[i]->box = 0;
      node[i]->color = 0;
    }
  
  if (bc->treeReadDistancesOn)
    {
      treeReadDistances(bc, pairmtx);
    }
  else
    {
      /* Calculate pairwise distance matrix */
      arrayOrder(bc->alignArr);
      
      for (i = 0; i < bc->alignArr->len - 1; ++i)
        {
          ALN *aln_i = &g_array_index(bc->alignArr, ALN, i);
          
          int j = i+1;
          for (j = i+1; j < bc->alignArr->len; ++j)
            {
              ALN *aln_j = &g_array_index(bc->alignArr, ALN, j);
              pairmtx[i][j] = 100.0 - identity(aln_i->seq, aln_j->seq, bc->penalize_gaps);
              
              if (bc->treeDistCorr == KIMURA) 
                pairmtx[i][j] = treeKimura(pairmtx[i][j]);
              else if (bc->treeDistCorr == JUKESCANTOR) 
                pairmtx[i][j] = treeJUKESCANTOR(pairmtx[i][j]);
              else if (bc->treeDistCorr == STORMSONN) 
                pairmtx[i][j] = treeSTORMSONN(pairmtx[i][j]);
              else if (bc->treeDistCorr == SCOREDIST) 
                pairmtx[i][j] = treeSCOREDIST(aln_i->seq, aln_j->seq, bc);
            }
        }
    }

  if (bc->treePrintDistances) 
    {
      double dist;
      
      for (i = 0; i < bc->alignArr->len; i++) 
        {
          if (!bc->treeCoordsOn) 
            {
              printf("%s\t", g_array_index(bc->alignArr, ALN, i).name);
            }
          else
            {
              printf("%s/%d-%d\t",
                     g_array_index(bc->alignArr, ALN, i).name,
                     g_array_index(bc->alignArr, ALN, i).start,
                     g_array_index(bc->alignArr, ALN, i).end);
            }
        }
      printf ("\n");
      
      for (i = 0; i < bc->alignArr->len; ++i) 
        {
          int j = 0;
          for (j = 0; j < bc->alignArr->len; ++j) 
            {
              if (i == j) 
                dist = 0.0;
              else if (i < j)
                dist = pairmtx[i][j];
              else
                dist = pairmtx[j][i];
              
              printf("%7.3f\t", dist);
            }
          printf ("\n");
        }
      exit(0);
    }
  
  /* Construct tree */
  int iter = 0;
  for (iter = 0; iter < bc->alignArr->len - 1; ++iter)
    {
      if (bc->treeMethod == NJ)
        {	
          int count = 0;
          
          /* Calculate vector r (avgdist) */
          for (i = 0; i < bc->alignArr->len; ++i)
            {
              if (!node[i]) 
                continue;
              
              avgdist[i] = 0.0;
              
              int j = 0;
              for (count=0, j = 0; j < bc->alignArr->len; j++)
                {
                  if (!node[j]) 
                    continue;
                  
                  avgdist[i] += (i < j ? pairmtx[i][j] : pairmtx[j][i]);
                  count++;
                }
              
              if (count == 2)	/* Hack, to accommodate last node */
                avgdist[i] = 1;
              else
                avgdist[i] /= 1.0*(count - 2);
            }
          
          /* Calculate corrected matrix Dij */
          if (1 /* gjm */)
            {
              for (i = 0; i < bc->alignArr->len - 1; ++i)
                {
                  if (!node[i])
                    continue;
                  
                  int j = i+1;
                  for (j = i+1; j < bc->alignArr->len; ++j)
                    {
                      if (!node[j]) 
                        continue;
                      
                      Dmtx[i][j] = pairmtx[i][j] - (avgdist[i] + avgdist[j]);
                    }
                }
            }
          else
            {		/* Felsenstein */
              double Q = 0.0;
              
              for (i = 0; i < bc->alignArr->len - 1; ++i)
                {
                  if (!node[i])
                    continue;
                  
                  int j = i+1;
                  for (j = i+1; j < bc->alignArr->len; j++)
                    {
                      if (!node[j]) 
                        continue;
                      
                      Q += pairmtx[i][j];
                    }
                }
              
              for (i = 0; i < bc->alignArr->len - 1; ++i)
                {
                  if (!node[i])
                    continue;
                  
                  int j = i+1;
                  for (j = i+1; j < bc->alignArr->len; j++)
                    {
                      if (!node[j])
                        continue;
                      
                      Dmtx[i][j] = (pairmtx[i][j] +
                                    (2.0*Q)/(count-(count == 2 ? 1 : 2)) - 
                                    avgdist[i] - avgdist[j]) / 2.0;
                    }
                }
            }
          curMtx = Dmtx;
        }
      else 
        {
          curMtx = pairmtx;
        }
      
#ifdef DEBUG
      {
        printf("Node status, Avg distance:\n");
        for (i = 0; i < bc->alignArr->len; i++)
          printf("%6d ", (node[i] ? 1 : 0)); 
        printf("\n");
        for (i = 0; i < bc->alignArr->len; i++)
          printf("%6.2f ", avgdist[i]); 
        printf("\n\nPairdistances, corrected pairdistances:");
        printMtx(pairmtx);
        printMtx(Dmtx);
        printf("\n");
      }
#endif
      
      /* Find smallest distance pair in pairmtx */
      maxi = -1;
      maxj = -1;
      maxid = 1000000;
      pmaxid = 1000000;

      for (i = 0; i < bc->alignArr->len - 1; ++i) 
        {
          if (!node[i]) 
            continue;
          
          int j = i+1;
          for (j = i+1; j < bc->alignArr->len; ++j) 
            {
              if (!node[j]) 
                continue;
              
              /* printf("iter %d, i=%d. j=%d, dist= %f\n", iter, i, j, curMtx[i][j]);*/
              if (curMtx[i][j] < maxid) 
                {
                  maxid = curMtx[i][j];
                  pmaxid = pairmtx[i][j];
                  maxi = i;
                  maxj = j;
                }
              else if (bc->treeMethod == NJ && Dmtx[i][j] == maxid && pairmtx[i][j] < pmaxid) 
                {
                  /* To resolve ties - important for tree look! */
                  maxi = i;
                  maxj = j;
                }
            }
        }
      
      maxid = pairmtx[maxi][maxj]; /* Don't want to point to Dmtx in NJ */

      /* Merge rows & columns of maxi and maxj into maxi
       Recalculate distances to other nodes */
      for (i = 0; i < bc->alignArr->len; ++i)
        {
          if (!node[i]) 
            continue;
          
          if (i < maxi) 
            trg = &pairmtx[i][maxi];
          else if (i > maxi) 
            trg = &pairmtx[maxi][i];
          else continue;
          
          if (i < maxj) 
            src = &pairmtx[i][maxj];
          else if (i > maxj) 
            src = &pairmtx[maxj][i];
          else continue;
          
          if (bc->treeMethod == UPGMA) 
            *trg = (*trg + *src) / 2.0;
          else
            *trg = (*trg + *src - maxid) / 2.0;
        }
      
      /* Create node for maxi and maxj */
      newnode = createEmptyTreeNode();
      
      if (bc->treeMethod == UPGMA)
        {
          /* subtract lower branch lengths from absolute distance
           Horribly ugly, only to be able to share code UPGMA and NJ */
          TreeNode *tmpnode = node[maxi]->left;
          
          for (llen = maxid; tmpnode; tmpnode = tmpnode->left)
            llen -= tmpnode->branchlen;
          
          tmpnode = node[maxj]->right;
          
          for (rlen = maxid; tmpnode; tmpnode = tmpnode->right)
            rlen -= tmpnode->branchlen;
        }
      else
        {
          llen = (maxid + avgdist[maxi] - avgdist[maxj]) / 2.0;
          rlen = maxid - llen;
          
          if (iter == bc->alignArr->len - 2)
            { /* Root node */
              
              /* Not necessary anymore, the tree is re-balanced at the end which calls this too
               treeBalanceByWeight(node[maxi], node[maxj], &llen, &rlen);*/
              
              /* Put entire length of root branch in one leg so the rebalancing 
               will work properly (otherwise it is hard to take this branch into account */
              rlen += llen;
              llen = 0;
            }
        }
      
#ifdef DEBUG
      {
        printf("Iter %d: Merging %d and %d, dist= %f\n", 
               iter, maxi+1, maxj+1, curMtx[maxi][maxj]);
        printf("maxid= %f  llen= %f  rlen= %f\n", maxid, llen, rlen);
        printf("avgdist[left]= %f  avgdist[right]= %f\n\n", 
               avgdist[maxi], avgdist[maxj]);
      }
#endif
      
      newnode->left = node[maxi];	
      newnode->left->branchlen = llen;
      
      newnode->right = node[maxj];
      newnode->right->branchlen = rlen;
      
      newnode->organism = (node[maxi]->organism == node[maxj]->organism ?
                           node[maxi]->organism : NULL);
      
      node[maxi] = newnode;
      node[maxj] = NULL;
    }

  fillParents(newnode, newnode->left);
  fillParents(newnode, newnode->right);
  
#ifdef DEBUG
  treeDraw(newnode) ;
#endif
  
  if (bc->treeMethod == UPGMA)
    newnode->branchlen = 100 - maxid ;
  
  if (bc->treeMethod == NJ)
    newnode = treeFindBalance(bc, newnode) ;
  
  fillOrganism(newnode);
  
  if (doBootstrap && bc->treebootstraps) 
    treeBootstrapStats(bc, newnode);
  
  return newnode ;
}


/***********************************************************
 *                        Drawing                          *
 ***********************************************************/

/* Clear any cached drawables and redraw everything */
static void redrawBelvuTree(GtkWidget *belvuTree)
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
  
  widgetClearCachedDrawable(properties->treeArea, NULL);
  
  gtk_widget_queue_draw(belvuTree);
}

/* Draws clickable boxes for horizontal lines.  The routine must be in sync with
 * treeDrawNode, but can not be integrated since a box may
 * accidentally overwrite some text.  Therefore all boxes must be 
 * drawn before any text or lines.
 */
//static double treeDrawNodeBox(BelvuContext *bc, Tree *tree, TreeNode *node, double x) 
//{
//  double y, yl, yr;
//  int box;
//  
//  if (!node) 
//    return 0.0;
//  
//  yl = treeDrawNodeBox(bc, tree, node->left, x + node->branchlen*treeScale);
//  yr = treeDrawNodeBox(bc, tree, node->right, x + node->branchlen*treeScale);
//  
//  if (yl) 
//    {
//      y = (yl + yr) / 2.0;
//    }
//  else 
//    {
//      y = bc->tree_y++;
//    }
//  
//  /* Make box around horizontal lines */
//  //  box = graphBoxStart();
//  //  graphLine(x + node->branchlen*treeScale, y, x, y);
//  //  oldlinew = graphLinewidth(0.0); graphColor(BG);
//  //  graphRectangle(x + node->branchlen*treeScale, y-0.5, x, y+0.5);
//  //  graphColor(BLACK); graphLinewidth(oldlinew);
//  //  graphBoxEnd();
//  
//  //    graphAssociate(assVoid(100+box), node);
//  
//  node->box = box;
//  tree->lastNodeBox = box;
//  
//  return y;
//}


/* The actual tree drawing routine.
 Note: must be in sync with treeDrawNodeBox, which draws clickable
 boxes first.
 */
static double treeDrawNode(BelvuContext *bc, 
                           GtkWidget *widget,
                           GdkDrawable *drawable, 
                           GdkGC *gc, 
                           BelvuTreeProperties *properties,
                           GdkColor *defaultColor, 
                           Tree *tree, 
                           TreeNode *node, 
                           double x) 
{
  double y, yl, yr;
  
  if (!node) 
    return 0.0;
  
  const int curX = x + (node->branchlen * (properties->treeScale * properties->charWidth));

  GdkGC *leftGc = gdk_gc_new(drawable);
  gdk_gc_copy(leftGc, gc);
  yl = treeDrawNode(bc, widget, drawable, leftGc, properties, defaultColor, tree, node->left, curX);
  
  gdk_gc_set_foreground(gc, defaultColor);
  yr = treeDrawNode(bc, widget, drawable, gc, properties, defaultColor, tree, node->right, curX);
  
  if (yl) 
    {
      /* internal node */
      y = (yl + yr) / 2.0;
      
      /* connect children */
      gdk_draw_line(drawable, gc, curX, yr, curX, y);
      gdk_gc_copy(gc, leftGc);
      gdk_draw_line(drawable, gc, curX, yl, curX, y);
      
      if (node->left->organism != node->right->organism)
        {
          /* Reset color */
          gdk_gc_set_foreground(gc, defaultColor);
        }
    }
  else 
    {
      /* Sequence name */
      int box ;
      
      y = bc->tree_y * properties->charHeight;
      bc->tree_y++;
      
      if (bc->treeColorsOn && node->organism) 
        {
          /* Get the color for this organism */
          ALN aln;
          initAln(&aln);
          
          aln.organism = node->organism;
          
          int ip = 0;
          if (arrayFind(bc->organismArr, &aln, &ip, (void*)organism_order)) 
            {
              GdkColor color;
              int colorNum = g_array_index(bc->organismArr, ALN, ip).color;
              convertColorNumToGdkColor(colorNum, &color);
              
              gdk_gc_set_foreground(gc, &color);
	    }
	}
      else
        {
          /* Just use the default color */
          gdk_gc_set_foreground(gc, defaultColor);
        }
      
      /* Make clickable box for sequence */
      GdkGC *gcTmp = gdk_gc_new(drawable);
      gdk_gc_set_foreground(gcTmp, defaultColor);
      drawText(widget, drawable, gcTmp, curX + properties->charWidth, y - properties->charHeight / 2, node->name);
      //      graphAssociate(assVoid(100+box), node);
      
      if (bc->highlightedAln && node->aln == bc->highlightedAln) 
        {
          //graphBoxDraw(box, WHITE, BLACK);
          tree->currentPickedBox = box;
	}
      else if (node->aln) 
        {
          //graphBoxDraw(box, BLACK, node->aln->color);
	}	    
      
      if (properties->showOrganism && node->organism) 
        {
          drawText(widget, drawable, gc, curX + (2 + strlen(node->name)) * properties->charWidth, y - properties->charHeight / 2, node->organism);
        }
      
        {
        int pos = curX + strlen(node->name);
          
        if (pos > bc->maxTreeWidth) 
          bc->maxTreeWidth = pos;
        }
    }
  
  /* Horizontal branches */
  gdk_draw_line(drawable, gc, curX, y, x, y);
  
  if (properties->showBranchLen && node->branchlen) 
    {
      char *tmpStr = blxprintf("%.1f", node->branchlen);
      double pos = x + (node->branchlen - strlen(tmpStr)) * properties->treeScale * properties->charWidth * 0.5;

      drawText(widget, drawable, gc, pos, y, tmpStr);
      
      g_free(tmpStr);
    }
  
  if (bc->treebootstraps && !node->name && node != bc->treeHead  && !bc->treebootstrapsDisplay) 
    {
      GdkColor *color = getGdkColor(BELCOLOR_TREE_BOOTSTRAP, bc->defaultColors, FALSE, FALSE);
      gdk_gc_set_foreground(gc, color);

      char *tmpStr = blxprintf("%.0f", node->boot);
      double pos = curX - strlen(tmpStr) - 0.5;
      
      if (pos < 0.0) 
        pos = 0;
      
      printf("%f  %f   \n", node->boot, pos);
      drawText(widget, drawable, gc, pos, y, tmpStr);

      g_free(tmpStr);
      gdk_gc_set_foreground(gc, defaultColor);
    }
  
  g_object_unref(leftGc);

  return y;
}


static Tree* createEmptyTree()
{
  Tree *result = g_malloc(sizeof *result);
  
  result->head = NULL;
  result->lastNodeBox = 0;
  result->currentPickedBox = 0;
  
  return result;
}


static void drawBelvuTree(GtkWidget *widget, GdkDrawable *drawable, BelvuTreeProperties *properties)
{
  BelvuContext *bc = properties->bc;
  
  
//  int i;
//  double oldlinew;
  
  Tree *treeStruct = createEmptyTree();
  treeStruct->head = properties->treeHead;
  
  //                           (bc->treeMethod == UPGMA ? 130 : 110) / fontwidth * screenWidth, 
  //                           (bc->alignArr->len + 7) / fontheight * screenHeight);
  //  graphRegister(PICK, treeboxPick);
  
                    
  GdkGC *gc = gdk_gc_new(drawable);
  GdkColor *defaultColor = getGdkColor(BELCOLOR_TREE_DEFAULT, bc->defaultColors, FALSE, FALSE);
  gdk_gc_set_line_attributes(gc, properties->lineWidth * properties->charWidth, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
  
  bc->maxTreeWidth = 0;
  bc->tree_y = 1;
  treeDrawNode(bc, widget, drawable, gc, properties, defaultColor, treeStruct, treeStruct->head, properties->treeRect.x);

  
  //graphTextBounds(bc->maxTreeWidth+2 + (treeShowOrganism ? 25 : 0), nseq+6);
  
  /* Draw scale */
//  if (bc->treeMethod == UPGMA) 
//    {
//      bc->tree_y += 2;
//      //graphLine(1*treeScale, tree_y, 101*treeScale, tree_y);
//      for (i=1; i<=101; i += 10) 
//        {
//          //graphLine(i*treeScale, tree_y, i*treeScale, tree_y+0.5);
//          
//          if (i < 101) 
//            {
//              //graphLine((i+5)*treeScale, tree_y, (i+5)*treeScale, tree_y+0.5);
//            }
//          
//          //graphText(messprintf("%d", i-1), (i-0.5)*treeScale, tree_y+1);
//        }
//    }
//  else
//    {
//      bc->tree_y += 3;
//      //graphText("0.1", 5*treeScale, tree_y-1);
//      //graphLine(1*treeScale, tree_y, 11*treeScale, tree_y);
//      //graphLine(1*treeScale, tree_y-0.5, 1*treeScale, tree_y+0.5);
//      //graphLine(11*treeScale, tree_y-0.5, 11*treeScale, tree_y+0.5);
//    }
//  //graphLinewidth(oldlinew);
//  
//  if (bc->treeMethod == NJ) 
//    {	
//      double lweight, rweight;
//      TreeNode *tree = treeStruct->head;
//      
//      lweight = treeSize3way(tree->left, tree);
//      rweight = treeSize3way(tree->right, tree);
//      
//      //      graphText((debug ? messprintf("Tree balance = %.1f (%.1f-%.1f)", 
//      //                                    fabsf(lweight - rweight), lweight, rweight) :
//      //                 messprintf("Tree balance = %.1f", fabsf(lweight - rweight))),
//      //                14, tree_y-0.5);
//    }
//  
  //    graphRedraw();
  
  g_object_unref(gc);
}


/* Expose function for the drawing area. The main belvuTree widget is passed as
 * the user data. */
static gboolean onExposeBelvuTree(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *belvuTree = GTK_WIDGET(data);
  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, properties->treeRect.x * 2 + properties->treeRect.width, 
                                          properties->treeRect.y * 2 + properties->treeRect.height);
          drawBelvuTree(widget, bitmap, properties);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
        }
      else
	{
	  g_warning("Failed to draw Belvu tree [%p] - could not create bitmap.\n", widget);
	}
    }
  
  return TRUE;
}


/***********************************************************
 *                    Settings dialog                      *
 ***********************************************************/

/* Callback to set the value of a boolean pointer (passed as the user data)
 * when the given toggle button has changed value. */
static gboolean onBoolChangedCallback(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean *value = (gboolean*)data;
  *value = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  return TRUE;
}

/* Callback to set the value of an integer (passed as the user data) when the
 * given text entry (which must contain a double or int as text) has changed value */
static gboolean onDoubleChangedCallback(GtkWidget *entry, const gint responseId, gpointer data)
{
  double *value = (double*)data;
  
  const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
  *value = g_strtod(text, NULL);
  
  return TRUE;
}

/* Utility to create a check button and place it in the given parent box. The
 * button will control the setting of the given boolean value, which will be 
 * updated when the button's parent dialog gets a response. */
static void createCheckButton(const char *mnemonic, gboolean *value, GtkTable *table, const int row, const int col)
{
  g_assert(value);
  int pad = 2;
  
  GtkWidget *button = gtk_check_button_new_with_mnemonic(mnemonic);
  gtk_table_attach(table, button, col, col + 2, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, pad, pad);
  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), *value);
  widgetSetCallbackData(button, onBoolChangedCallback, value);
}

/* Utility to create a text entry box for updating the given integer value.
 * The text entry is added to the given box, which should belong to a dialog.
 * The value will be updated when the dialog gets a response. */
static void createDoubleTextEntry(const char *labelText, double *value, GtkTable *table, int row, int col)
{
  g_assert(value);
  int pad = 2;

  /* Create the label */
  GtkWidget *label = gtk_label_new(labelText);
  gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
  gtk_table_attach(table, label, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, pad, pad);
  
  /* Create the text entry */
  GtkWidget *entry = gtk_entry_new();
  gtk_table_attach(table, entry, col + 1, col + 2, row, row + 1, GTK_EXPAND | GTK_FILL, GTK_SHRINK, pad, pad);
  widgetSetCallbackData(entry, onDoubleChangedCallback, value);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

  char *defaultInput = convertDoubleToString(*value, 2);
  gtk_entry_set_text(GTK_ENTRY(entry), defaultInput);

  
  const int defaultLen = min(strlen(defaultInput) * 8, 500);
  gtk_widget_set_size_request(entry, defaultLen, -1);

  g_free(defaultInput);
}

static void createTreeDisplayOptsButtons(GtkBox *box, 
                                         double *treeScale,
                                         double *lineWidth,
                                         gboolean *showBranchLen, 
                                         gboolean *showOrganism)
{
  /* Put the display options in a table in a frame */
  GtkContainer *frame = GTK_CONTAINER(gtk_frame_new("Display options"));
  gtk_box_pack_start(box, GTK_WIDGET(frame), FALSE, FALSE, 12);
  
  GtkTable *table = GTK_TABLE(gtk_table_new(2, 2, FALSE));
  gtk_container_add(frame, GTK_WIDGET(table));

  createDoubleTextEntry("Tree scale: ", treeScale, table, 0, 0);
  createDoubleTextEntry("Line width: ", lineWidth, table, 1, 0);
  createCheckButton("Display branch lengths", showBranchLen, table, 2, 0);
  createCheckButton("Display organism", showOrganism, table, 3, 0);
}

static void createTreeInteractionButtons(GtkBox *box)
{
  
}

/* Utility function to create the content for the tree settings dialog */
void createTreeSettingsDialogContent(BelvuContext *bc, 
                                     GtkWidget *dialog, 
                                     double *treeScale,
                                     double *lineWidth,
                                     gboolean *showBranchLen, 
                                     gboolean *showOrganism)
{
  GtkBox *vbox = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), GTK_WIDGET(vbox), FALSE, FALSE, 0);
  
  createTreeDisplayOptsButtons(vbox, treeScale, lineWidth, showBranchLen, showOrganism);
  createTreeInteractionButtons(vbox);
}


/* Called when the user makes a response on the tree settings dialog. The tree
 * window is passed as the user data. */
void onResponseTreeSettingsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  GtkWidget *belvuTree = GTK_WIDGET(data);
  gboolean destroy = TRUE;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      /* Call all of the callbacks for each individual widget to update the 
       * properties. Then refresh the window. Destroy if successful. */
      destroy = widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      redrawBelvuTree(belvuTree);
      break;
      
    case GTK_RESPONSE_APPLY:
      /* Never destroy */
      destroy = FALSE;
      widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      redrawBelvuTree(belvuTree);
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      destroy = TRUE;
      break;
      
    default:
      break;
  };
  
  if (destroy)
    {
      gtk_widget_destroy(GTK_WIDGET(dialog));
    }
}


/* Dialog to allow the user to edit the settings for a tree */
void showTreeSettingsDialog(GtkWidget *belvuTree)
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
  BelvuContext *bc = properties->bc;
  
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Belvu - Tree Settings", 
                                       GTK_WINDOW(belvuTree), 
                                       GTK_DIALOG_DESTROY_WITH_PARENT | GTK_DIALOG_MODAL,
                                       GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                       GTK_STOCK_APPLY, GTK_RESPONSE_APPLY,
                                       GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                       NULL);
      
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseTreeSettingsDialog), belvuTree);
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_APPLY);
  
  createTreeSettingsDialogContent(bc, dialog, 
                                  &properties->treeScale, &properties->lineWidth,
                                  &properties->showBranchLen, &properties->showOrganism);
  
  gtk_widget_show_all(dialog);
  gtk_window_present(GTK_WINDOW(dialog));
}

/***********************************************************
 *                          Sizing                         *
 ***********************************************************/

static void calculateBelvuTreeBorders(GtkWidget *belvuTree)
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
  
  properties->treeRect.x = DEFAULT_XPAD;
  properties->treeRect.y = DEFAULT_YPAD;
  properties->treeRect.width = 1000; /* to do: calculate real width and height */
  properties->treeRect.height = 800;
  
  gtk_layout_set_size(GTK_LAYOUT(properties->treeArea), properties->treeRect.x * 2 + properties->treeRect.width, 
                      properties->treeRect.y * 2 + properties->treeRect.height);
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* Set style properties for the belvu alignment widgets */
static void setBelvuTreeStyle(BelvuContext *bc, GtkWidget *belvuTree)
{
  GdkColor *bgColor = getGdkColor(BELCOLOR_BACKGROUND, bc->defaultColors, FALSE, FALSE);
  
  gtk_widget_modify_bg(belvuTree, GTK_STATE_NORMAL, bgColor);
}


/* Create a widget for drawing a tree. Returns the GtkLayout widget that 
 * does the drawing. Places the widget into the given parent box. */
static GtkWidget* createBelvuTreeWidget(BelvuContext *bc, TreeNode *treeHead, GtkBox *box)
{
  /* Create the drawing area */
  GtkWidget *treeArea = gtk_layout_new(NULL, NULL);
  
  /* Make it scrollable. The scroll-window will be the outermost container and
   * therefore will be the widget that we treat as the "belvu tree". */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_box_pack_start(box, scrollWin, TRUE, TRUE, 0);
  
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  gtk_container_add(GTK_CONTAINER(scrollWin), treeArea);

  /* Set the style */
  setBelvuTreeStyle(bc, treeArea);
  
  return treeArea;
}


static void setTreeWindowStyleProperties(GtkWidget *window)
{
  gtk_widget_set_name(window, BELVU_TREE_WINDOW_NAME);
  
  /* Set the initial window size based on some fraction of the screen size */
  GdkScreen *screen = gtk_widget_get_screen(window);
  const int screenWidth = gdk_screen_get_width(screen);
  const int screenHeight = gdk_screen_get_height(screen);
  
  const int width = screenWidth * DEFAULT_TREE_WINDOW_WIDTH_FRACTION;
  const int height = screenHeight * DEFAULT_TREE_WINDOW_HEIGHT_FRACTION;
  gtk_window_set_default_size(GTK_WINDOW(window), width, height);
  
  /* Set the initial position */
  const int x = (screenWidth - width) / 4;
  const int y = (screenHeight - height) / 4;
  gtk_window_move(GTK_WINDOW(window), x, y);
}


/* Display an existing tree. Opens in a new window. */
void createBelvuTreeWindow(BelvuContext *bc, TreeNode *treeHead)
{
  /* Create the window */
  GtkWidget *belvuTree = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setTreeWindowStyleProperties(belvuTree);

  char *title = blxprintf("Belvu - %s using %s distances of %s", 
                          bc->treeMethod == NJ ? "Neighbor-joining tree" : "UPGMA tree",
                          bc->treeDistString,
                          bc->Title);
  
  gtk_window_set_title(GTK_WINDOW(belvuTree), title);
  g_free(title);
  
  /* Create the context menu and set a callback to show it */
  GtkUIManager *uiManager = createUiManager(belvuTree, bc, NULL);
  GtkWidget *contextmenu = createBelvuMenu(belvuTree, "/TreeContextMenu", uiManager);
  
  gtk_widget_add_events(belvuTree, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(belvuTree), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  
  /* We'll place everything in a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(belvuTree), vbox);
  
  /* Add the alignment section */
  GtkWidget *treeArea = createBelvuTreeWidget(bc, treeHead, GTK_BOX(vbox));

  /* Set the properties on the tree window */
  belvuTreeCreateProperties(belvuTree, bc, treeArea, treeHead);

  /* Set the initial size (must be called after properties are set) */
  calculateBelvuTreeBorders(belvuTree);

  /* Connect expose signal for the drawing area. Pass the window as the data. */
  g_signal_connect(G_OBJECT(treeArea), "expose-event", G_CALLBACK(onExposeBelvuTree), belvuTree);  

  gtk_widget_show_all(belvuTree);
  gtk_window_present(GTK_WINDOW(belvuTree));
}


/* Create a new tree and create a tree-window to display it. */
void createAndShowBelvuTree(BelvuContext *bc)
{
  separateMarkupLines(bc);
  
  bc->treeHead = treeMake(bc, TRUE);
  createBelvuTreeWindow(bc, bc->treeHead);
  
  reInsertMarkupLines(bc);
  
}

