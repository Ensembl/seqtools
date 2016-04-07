/*  File: blxcontext.cpp
 *  Author: Gemma Barson, 2016-04-06
 *  Copyright (c) 2016 Genome Research Ltd
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
 * Description: See blxcontext.hpp
 *----------------------------------------------------------------------------
 */


#include <algorithm>

#include <blixemApp/blxcontext.hpp>
#include <blixemApp/blixem_.hpp>
#include <seqtoolsUtils/blxGff3Parser.hpp>


using namespace std;



/***********************************************************
 *                     Internal functions                  *
 ***********************************************************/

namespace // unnamed namespace
{

/* This returns the name for the given flag enum. It returns null if the name
 * has not been set. This is used to save settings in the config file; only 
 * flags whose name is set in this function will be saved. */
static const char* getFlagName(const BlxFlag flag)
{
  /* Create an array of names for all of the relevant settings */
  static const char* names[BLXFLAG_NUM_FLAGS] = {0, 0};
  
  if (!names[1]) /* only populate it once */
    {
      names[BLXFLAG_HIGHLIGHT_DIFFS] = SETTING_NAME_HIGHLIGHT_DIFFS;
      names[BLXFLAG_HIGHLIGHT_VARIATIONS] = SETTING_NAME_HIGHLIGHT_VARIATIONS;
      names[BLXFLAG_SHOW_VARIATION_TRACK] = SETTING_NAME_SHOW_VARIATION_TRACK;
      names[BLXFLAG_SHOW_UNALIGNED] = SETTING_NAME_SHOW_UNALIGNED;
      names[BLXFLAG_SHOW_UNALIGNED_SELECTED] = SETTING_NAME_SHOW_UNALIGNED_SELECTED;
      names[BLXFLAG_LIMIT_UNALIGNED_BASES] = SETTING_NAME_LIMIT_UNALIGNED_BASES;
      names[BLXFLAG_SHOW_POLYA_SITE] = SETTING_NAME_SHOW_POLYA_SITE;
      names[BLXFLAG_SHOW_POLYA_SITE_SELECTED] = SETTING_NAME_SHOW_POLYA_SITE_SELECTED;
      names[BLXFLAG_SHOW_POLYA_SIG] = SETTING_NAME_SHOW_POLYA_SIG;
      names[BLXFLAG_SHOW_POLYA_SIG_SELECTED] = SETTING_NAME_SHOW_POLYA_SIG_SELECTED;
      names[BLXFLAG_SHOW_SPLICE_SITES] = SETTING_NAME_SHOW_SPLICE_SITES;
      names[BLXFLAG_SHOW_MAYBE_CANONICAL] = SETTING_NAME_SHOW_MAYBE_CANONICAL;
      names[BLXFLAG_SHOW_COLINEARITY] = SETTING_NAME_SHOW_COLINEARITY;
      names[BLXFLAG_SHOW_COLINEARITY_SELECTED] = SETTING_NAME_SHOW_COLINEARITY_SELECTED;
    }

  if (flag <= BLXFLAG_MIN || flag >= BLXFLAG_NUM_FLAGS)
    return NULL;
  else
    return names[flag];
}


/* Whether to include the given msp type in depth coverage calculations */
static gboolean includeTypeInCoverage(BlxMspType mspType)
{
  return (mspType == BLXMSP_MATCH);
}


/* Utility to get the depth-counter enum from the given character */
DepthCounter getDepthCounterForChar(const char c, const BlxStrand strand)
{
  DepthCounter result = DEPTHCOUNTER_NONE;

  switch (c)
    {
    case 'a': //fall through
    case 'A': 
      if (strand == BLXSTRAND_REVERSE)
        result = DEPTHCOUNTER_A_R;
      else
        result = DEPTHCOUNTER_A_F;
      break;
    case 'c': //fall through 
    case 'C': 
      if (strand == BLXSTRAND_REVERSE)
        result = DEPTHCOUNTER_C_R;
      else
        result = DEPTHCOUNTER_C_F;
      break;
    case 'g':  //fall through
    case 'G': 
      if (strand == BLXSTRAND_REVERSE)
        result = DEPTHCOUNTER_G_R; 
      else
        result = DEPTHCOUNTER_G_F; 
      break;
    case 'u':  //fall through
    case 'U':   //fall through
    case 't':  //fall through
    case 'T': 
      if (strand == BLXSTRAND_REVERSE)
        result = DEPTHCOUNTER_T_R; 
      else
        result = DEPTHCOUNTER_T_F; 
      break;
    case 'n':  //fall through
    case 'N': 
      if (strand == BLXSTRAND_REVERSE)
        result = DEPTHCOUNTER_N_R; 
      else
        result = DEPTHCOUNTER_N_F; 
      break;
    case '.': // indicates a gap
      if (strand == BLXSTRAND_REVERSE)
        result = DEPTHCOUNTER_GAP_R; 
      else
        result = DEPTHCOUNTER_GAP_F; 
      break;
    default:
      break;
    }

  return result;
}


/* utility to free the given pointer and set it to null */
static void freeAndNull(gpointer *ptr)
{
  if (ptr && *ptr)
    {
      g_free(*ptr);
      *ptr = NULL;
    }
}


} // unnamed namespace


/***********************************************************
 *                     BlxContext class                    *
 ***********************************************************/

BlxContext::BlxContext(CommandLineOptions *options,
                       const IntRange* const refSeqRange_in,
                       const IntRange* const fullDisplayRange_in,
                       const char *paddingSeq_in,
                       GArray* featureLists_in[],
                       GList *seqList_in,
                       GSList *supportedTypes_in,
                       GtkWidget *widget_in,
                       GtkWidget *statusBar_in,
                       const gboolean External_in,
                       GSList *styles_in)
{
  statusBar = statusBar_in;
  
  refSeq = options->refSeq;
  refSeqName = options->refSeqName ? g_strdup(options->refSeqName) : g_strdup("Blixem-seq");
  refSeqRange.set(refSeqRange_in);
  fullDisplayRange.set(fullDisplayRange_in);
  refSeqOffset = options->refSeqOffset;
  optionalColumns = options->optionalColumns;

  mspList = options->mspList;
  columnList = options->columnList;
  styles = styles_in;
  
  int typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    {
      featureLists[typeId] = featureLists_in[typeId];
    }
  
  geneticCode = options->geneticCode;
  blastMode = options->blastMode;
  seqType = options->seqType;
  numFrames = options->numFrames;
  paddingSeq = paddingSeq_in;
  bulkFetchDefault = options->bulkFetchDefault;
  userFetchDefault = options->userFetchDefault;
  optionalFetchDefault = options->optionalFetchDefault;
  fetchMethods = options->fetchMethods;
  dataset = g_strdup(options->dataset);
  matchSeqs = seqList_in;
  supportedTypes = supportedTypes_in;
  
  displayRev = FALSE;
  external = External_in;
  
  selectedSeqs = NULL;
  sequenceGroups = NULL;
  matchSetGroup = NULL;
  
  dotterRefType = BLXDOTTER_REF_AUTO;
  dotterMatchType = BLXDOTTER_MATCH_SELECTED;
  dotterAdhocSeq = NULL;
  dotterHsps = FALSE;
  dotterSleep = FALSE;
  dotterStart = UNSET_INT;
  dotterEnd = UNSET_INT;
  dotterZoom = 0;
  
  defaultColors = NULL;
  usePrintColors = FALSE;
  windowColor = options->windowColor;
  
  createColors(widget_in);
  
  initialiseFlags(options);
    
  /* Null out all the entries in the dialogs list */
  int dialogId = 0;
  for ( ; dialogId < BLXDIALOG_NUM_DIALOGS; ++dialogId)
    {
      dialogList[dialogId] = NULL;
    }
    
  spawnedProcesses = NULL;
  minDepth = 0;
  maxDepth = 0;

  for (int counter = DEPTHCOUNTER_NONE + 1; counter < DEPTHCOUNTER_NUM_ITEMS; ++counter)
    depthArray[counter] = NULL;
 
  loadSettings();

  /* do this after loading settings because the passed-in squashed 
   * matches option should override the saved option in the settings */
  modelId = options->squashMatches ? BLXMODEL_SQUASHED : BLXMODEL_NORMAL;

  fetch_debug = options->fetch_debug;

#ifdef PFETCH_HTML
  ipresolve = options->ipresolve;
  cainfo = options->cainfo;
#endif

  /* Calculate the font size */
  if (widget_in)
    getFontCharSize(widget_in, widget_in->style->font_desc, &m_charWidth, &m_charHeight);

}


BlxContext::~BlxContext()
{
  /* Free allocated strings */
  freeAndNull((gpointer*)(&dataset));
  freeAndNull((gpointer*)(&refSeqName));

  /* Free table of fetch methods and the fetch-method structs */
  /* to do */
      
  /* Free the list of selected sequence names (not the names themselves
   * because we don't own them). */
  if (selectedSeqs)
    {
      g_list_free(selectedSeqs);
      selectedSeqs = NULL;
    }
      
  deleteAllSequenceGroups();

  /* Free the color array */
  if (defaultColors)
    {
      int i = BLXCOLOR_MIN + 1;
      for (; i < BLXCOL_NUM_COLORS; ++i)
        {
          BlxColor *blxColor = &g_array_index(defaultColors, BlxColor, i);
          destroyBlxColor(blxColor);
        }

      g_array_free(defaultColors, TRUE);
      defaultColors = NULL;
    }

  /* destroy the feature lists. note that the stored msps are owned
   * by the msplist, not by the feature lists */
  int typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    g_array_free(featureLists[typeId], FALSE);
      
  destroyMspList(&(mspList));
  destroyBlxSequenceList(&(matchSeqs));
  blxDestroyGffTypeList(&(supportedTypes));
  killAllSpawned();
}


/* Free the memory used by the given sequence group and its members. */
void BlxContext::destroySequenceGroup(SequenceGroup **seqGroup)
{
  if (seqGroup && *seqGroup)
    {
      /* Remove it from the list of groups */
      sequenceGroups = g_list_remove(sequenceGroups, *seqGroup);
      
      /* If this is pointed to by the match-set pointer, null it */
      if (*seqGroup == matchSetGroup)
        {
          matchSetGroup = NULL;
        }
      
      /* Free the memory used by the group name */
      if ((*seqGroup)->groupName)
        {
          g_free((*seqGroup)->groupName);
        }
      
      /* Free the list of sequences */
      if ((*seqGroup)->seqList)
        {
          freeStringList(&(*seqGroup)->seqList, (*seqGroup)->ownsSeqNames);
        }
      
      g_free(*seqGroup);
      *seqGroup = NULL;
    }
}


/* Delete all groups */
void BlxContext::deleteAllSequenceGroups()
{
  GList *groupItem = sequenceGroups;
  
  for ( ; groupItem; groupItem = groupItem->next)
    {
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      destroySequenceGroup(&group);
    }
  
  g_list_free(sequenceGroups);
  sequenceGroups = NULL;
  
  /* Reset the hide-not-in-group flags otherwise we'll hide everything! */
  flags[BLXFLAG_HIDE_UNGROUPED_SEQS] = FALSE ;
  flags[BLXFLAG_HIDE_UNGROUPED_FEATURES] = FALSE ;
}


double BlxContext::charWidth() const
{
  return m_charWidth;
}

double BlxContext::charHeight() const
{
  return m_charHeight;
}


BlxStrand BlxContext::activeStrand() const
{
  BlxStrand result = BLXSTRAND_NONE;

  result = displayRev ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD;
  
  return result;
}


/* Create the colors that blixem will use for various specific purposes */
void BlxContext::createColors(GtkWidget *widget)
{
  /* Initialise the array with empty BlxColor structs */
  defaultColors = g_array_sized_new(FALSE, FALSE, sizeof(BlxColor), BLXCOL_NUM_COLORS);
  int i = BLXCOLOR_MIN + 1;
  
  for ( ; i < BLXCOL_NUM_COLORS; ++i)
    {
      BlxColor *blxColor = (BlxColor*)g_malloc(sizeof(BlxColor));
      blxColor->name = NULL;
      blxColor->desc = NULL;
      g_array_append_val(defaultColors, *blxColor);
    }
  
  /* Get the default background color of our widgets (i.e. that inherited from the theme).
   * Convert it to a string so we can use the same creation function as the other colors */
  char *defaultBgColorStr = convertColorToString(&widget->style->bg[GTK_STATE_NORMAL]);
  createBlxColor(defaultColors, BLXCOLOR_BACKGROUND, "Background", "Background color", defaultBgColorStr, BLX_WHITE, "#bdbdbd", NULL);
  
  /* reference sequence */
  createBlxColor(defaultColors, BLXCOLOR_REF_SEQ, "Reference sequence", "Default background color for the reference sequence", BLX_YELLOW, BLX_VERY_LIGHT_GREY, BLX_DARK_YELLOW, NULL);
  
  /* matches */
  createBlxColor(defaultColors, BLXCOLOR_MATCH, "Exact match", "Exact match", BLX_LIGHT_CYAN, BLX_LIGHT_GREY, BLX_CYAN, NULL);
  createBlxColor(defaultColors, BLXCOLOR_CONS, "Conserved match", "Conserved match", BLX_VIOLET, BLX_VERY_LIGHT_GREY, BLX_DARK_VIOLET, NULL);
  createBlxColor(defaultColors, BLXCOLOR_MISMATCH, "Mismatch", "Mismatch", "#FFFFFF", BLX_WHITE, "#FED4EA", NULL);
  createBlxColor(defaultColors, BLXCOLOR_INSERTION, "Insertion", "Insertion", "#9E00FF", BLX_VERY_DARK_GREY, NULL, NULL);
  
  /* exons */
  createBlxColor(defaultColors, BLXCOLOR_EXON_START, "Exon start", "Exon start boundary", BLX_BLUE, BLX_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_EXON_END, "Exon end", "Exon end boundary", BLX_DARK_BLUE, BLX_GREY, NULL, NULL);

  createBlxColor(defaultColors, BLXCOLOR_EXON_FILL, "Exon fill color", "Exon fill color in big picture", BLX_PALE_YELLOW, BLX_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_EXON_LINE, "Exon line color", "Exon line color in big picture", BLX_BLUE, BLX_VERY_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_CDS_FILL, "CDS fill color", "Coding section fill color in big picture", BLX_LIGHT_GREEN, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_CDS_LINE, "CDS line color", "Coding section line color in big picture", BLX_DARK_GREEN, BLX_DARK_GREY, BLX_VERY_DARK_GREEN, NULL);
  createBlxColor(defaultColors, BLXCOLOR_UTR_FILL, "Exon fill color (UTR)", "Untranslated region fill color in big picture", BLX_LIGHT_RED, BLX_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_UTR_LINE, "Exon line color (UTR)", "Untranslated region line color in big picture", BLX_DARK_RED, BLX_VERY_DARK_GREY, BLX_VERY_DARK_RED, NULL);
  createBlxColor(defaultColors, BLXCOLOR_PARTIAL_EXON_CROSSHATCH, "Cross-hatch line color for partial exons", "Line color of cross-hatch highlighting for partial exons", BLX_GREY, BLX_GREY, NULL, NULL);
  
  /* codons */
  createBlxColor(defaultColors, BLXCOLOR_CODON, "Codon nucleotides", "Codon nucleotides", BLX_SKY_BLUE, BLX_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_MET, "MET codons", "MET codons", BLX_LAWN_GREEN, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_STOP, "STOP codons", "MET codons", BLX_SALMON_RED, BLX_LIGHT_GREY, NULL, NULL);
  
  /* SNPs */
  createBlxColor(defaultColors, BLXCOLOR_SNP, "SNPs", "SNPs", BLX_ORANGE, BLX_GREY, NULL, NULL);

  /* Big Picture */
  createBlxColor(defaultColors, BLXCOLOR_GRID_LINE, "Grid lines", "Big Picture grid lines", BLX_YELLOW, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_GRID_TEXT, "Grid text", "Big Picture grid text", BLX_BLACK, BLX_BLACK, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_HIGHLIGHT_BOX, "Highlight box", "Highlight box in the big picture", BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_PREVIEW_BOX, "Preview box", "Preview box in the big picture", BLX_BLACK, BLX_BLACK, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_MSP_LINE, "Big picture match line", "Color of the lines representing matches in the Big Picture", BLX_BLACK, BLX_BLACK, BLX_CYAN, BLX_GREY);

  /* groups */
  createBlxColor(defaultColors, BLXCOLOR_GROUP, "Default group color", "Default highlight color for a new group", BLX_ORANGE_RED, BLX_VERY_LIGHT_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_MATCH_SET, "Default match set color", "Default color for the match set group (applies only when it is created for the first time or after being deleted)", BLX_RED, BLX_VERY_LIGHT_GREY, NULL, NULL);

  /* colinearity */
  createBlxColor(defaultColors, BLXCOLOR_COLINEAR_PERFECT, "Perfect colinearity", "Color of lines joining alignment blocks with perfect colinearity", BLX_DARK_GREEN, BLX_LIGHT_GREY, BLX_DARK_GREEN, BLX_LIGHT_GREY);
  createBlxColor(defaultColors, BLXCOLOR_COLINEAR_IMPERFECT, "Imperfect colinearity", "Color of lines joining alignment blocks with imperfect colinearity", BLX_ORANGE, BLX_GREY, BLX_ORANGE, BLX_GREY);
  createBlxColor(defaultColors, BLXCOLOR_COLINEAR_NOT, "Not colinear", "Color of lines joining alignment blocks that are not colinear", BLX_RED, BLX_DARK_GREY, BLX_RED, BLX_DARK_GREY);

  /* polyA features */
  createBlxColor(defaultColors, BLXCOLOR_POLYA_TAIL, "polyA tail", "polyA tail", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_POLYA_SIGNAL, "polyA signal", "polyA signal", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_POLYA_SIGNAL_ANN, "Annotated polyA signal", "Annotated polyA signal", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_POLYA_SITE_ANN, "Annotated polyA site", "Annotated polyA site", BLX_RED, BLX_DARK_GREY, NULL, NULL);

  /* misc */
  createBlxColor(defaultColors, BLXCOLOR_UNALIGNED_SEQ, "Unaligned sequence", "Addition sequence in the match that is not part of the alignment", "#FFC432", BLX_WHITE, "#FFE8AD", NULL);
  createBlxColor(defaultColors, BLXCOLOR_CANONICAL, "Canonical intron bases", "The two bases at the start/end of the intron for the selected MSP are colored this color if they are canonical", BLX_GREEN, BLX_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_NON_CANONICAL, "Non-canonical intron bases", "The two bases at the start/end of the intron for the selected MSP are colored this color if they are not canonical", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_MAYBE_CANONICAL, "\"Maybe\" canonical intron bases", "The two bases at the start/end of the intron for the selected MSP are colored this color if they would be canonical if they were on the other strand", BLX_ORANGE, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_TREE_GRID_LINES, "Tree grid lines", "Tree grid lines", BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY);
  createBlxColor(defaultColors, BLXCOLOR_CLIP_MARKER, "Clipped-match indicator", "Marker to indicate a match has been clipped to the display range", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_COVERAGE_PLOT, "Coverage plot", "Coverage plot", BLX_ROYAL_BLUE, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_ASSEMBLY_GAP, "Assembly gaps", "Highlight color for assembly gaps", "#D14553", BLX_DARK_GREY, NULL, NULL);
  createBlxColor(defaultColors, BLXCOLOR_SELECTION, "Selection color", "Highlight color for selections", BLX_DARK_GREY, BLX_DARK_GREY, NULL, NULL);
  
  g_free(defaultBgColorStr);
}
/* Called on startup to set the initial state of the flags. Gets the state for
 * the settings from the config file if specified, otherwises uses hard-coded
 * defaults. */
void BlxContext::initialiseFlags(CommandLineOptions *options)
{
  /* Initialise all the flags to false */
  int flag = BLXFLAG_MIN + 1;
  for ( ; flag < BLXFLAG_NUM_FLAGS; ++flag)
    {
      flags[flag] = FALSE;
    }
  
  /* Set any specific flags that we want initialised to TRUE */
  flags[BLXFLAG_LIMIT_UNALIGNED_BASES] = TRUE;
  flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED] = TRUE;
  flags[BLXFLAG_SHOW_POLYA_SIG_SELECTED] = TRUE;
  flags[BLXFLAG_SHOW_SPLICE_SITES] = TRUE;
  flags[BLXFLAG_NEGATE_COORDS] = options->negateCoords;
  flags[BLXFLAG_HIGHLIGHT_DIFFS] = options->highlightDiffs;
  flags[BLXFLAG_SAVE_TEMP_FILES] = options->saveTempFiles;
  flags[BLXFLAG_ABBREV_TITLE] = options->abbrevTitle;
}


/* load settings from the config file */
void BlxContext::loadSettings()
{
  /* Override the defaults settings with those given in the config file, if any */
  GKeyFile *key_file = blxGetConfig();

  if (!key_file)
    return;
  
  GError *error = NULL;

  /* squash-matches */
  int squashMatches = g_key_file_get_integer(blxGetConfig(), SETTINGS_GROUP, SETTING_NAME_SQUASH_MATCHES, &error);
  
  if (error)
    {
      /* we don't care if it wasn't found; just clear the error */
      g_error_free(error);
      error = NULL;
    }
  else
    {
      modelId = squashMatches ? BLXMODEL_SQUASHED : BLXMODEL_NORMAL;
    }

  /* loop through all the flags and see if any of them are given */
  int flag = BLXFLAG_MIN + 1;
  
  for ( ; flag < BLXFLAG_NUM_FLAGS; ++flag)
    {
      const char *flagName = getFlagName((BlxFlag)flag);
      
      if (flagName)
        {
          int result = g_key_file_get_integer(key_file, SETTINGS_GROUP, flagName, &error);
          
          if (error)
            {
              g_error_free(error);
              error = NULL;
            }
          else
            {
              flags[flag] = result;
            }
        }
    }
}


/* Called by saveBlixemSettings; does the work to save the boolean flags */
void BlxContext::saveSettingsFlags(GKeyFile *key_file)
{
  /* loop through each (save-able) setting */
  int flag = BLXFLAG_MIN + 1;
  
  for ( ; flag < BLXFLAG_NUM_FLAGS; ++flag)
    {
      const char *flagName = getFlagName((BlxFlag)flag);

      if (flagName)
        g_key_file_set_integer(key_file, SETTINGS_GROUP, flagName, flags[flag]);
    }
}


/* Kill all processes spawned from this blixem */
void BlxContext::killAllSpawned()
{
  GSList *processes = spawnedProcesses;
  
  for ( ; processes; processes = processes->next)
    {
      pid_t pid = GPOINTER_TO_INT(processes->data);
      kill(pid, 9);
    }
    
  if (spawnedProcesses)
    {
      g_slist_free(spawnedProcesses);
      spawnedProcesses = NULL;
    }
}


/* Calculate the depth of coverage of short-reads for each reference sequence display coord.
 * depthArray must be the same length as displayRange. */
void BlxContext::calculateDepth(const int numUnalignedBases)
{
  /* Allocate the depth array, if null */
  const int displayLen = fullDisplayRange.length();
  
  if (displayLen < 1)
    return; 
  
  depthArray[DEPTHCOUNTER_ALL_F] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_GAP_F] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_A_F] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_C_F] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_G_F] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_T_F] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_N_F] = (int*)g_malloc0(sizeof(int) * displayLen);
  
  depthArray[DEPTHCOUNTER_ALL_R] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_GAP_R] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_A_R] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_C_R] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_G_R] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_T_R] = (int*)g_malloc0(sizeof(int) * displayLen);
  depthArray[DEPTHCOUNTER_N_R] = (int*)g_malloc0(sizeof(int) * displayLen);
  
  /* Initialise each entry to zero */  
  int i = 0;
  for ( ; i < displayLen; ++i)
    {
      depthArray[DEPTHCOUNTER_ALL_F][i] = 0;
      depthArray[DEPTHCOUNTER_GAP_F][i] = 0;
      depthArray[DEPTHCOUNTER_A_F][i] = 0;
      depthArray[DEPTHCOUNTER_C_F][i] = 0;
      depthArray[DEPTHCOUNTER_G_F][i] = 0;
      depthArray[DEPTHCOUNTER_T_F][i] = 0;
      depthArray[DEPTHCOUNTER_N_F][i] = 0;

      depthArray[DEPTHCOUNTER_ALL_R][i] = 0;
      depthArray[DEPTHCOUNTER_GAP_R][i] = 0;
      depthArray[DEPTHCOUNTER_A_R][i] = 0;
      depthArray[DEPTHCOUNTER_C_R][i] = 0;
      depthArray[DEPTHCOUNTER_G_R][i] = 0;
      depthArray[DEPTHCOUNTER_T_R][i] = 0;
      depthArray[DEPTHCOUNTER_N_R][i] = 0;
    }
  
  /* Loop through all MSP lists */
  int mspType = 0;
  
  for ( ; mspType < BLXMSP_NUM_TYPES; ++mspType)
    {
      /* Only include MSPs of relevant types */
      if (!includeTypeInCoverage((BlxMspType)mspType))
        continue;
      
      /* Loop through all MSPs in this list */
      GArray *mspArray = featureLists[mspType];
      const int fullDisplayLen = fullDisplayRange.length();
    
      i = 0;
      const MSP *msp = mspArrayIdx(mspArray, i);
  
      for ( ; msp; msp = mspArrayIdx(mspArray, ++i))
        {
          /* For each ref-seq coord that this alignment spans, increment the depth */
          int alignIdx = msp->displayRange.min();
          int qIdx = msp->qRange.min();

          for ( ; alignIdx <= msp->displayRange.max(); ++alignIdx, ++qIdx)
            {
              /* Convert the msp coord to a zero-based coord. Note that parts of the
               * msp range may be outside the ref seq range. */
              const int displayIdx = alignIdx - fullDisplayRange.min();

              if (displayIdx >= 0 && displayIdx < fullDisplayLen)
                {
                  /* Increment the main counter */
                  if (msp->qStrand == BLXSTRAND_REVERSE)
                    depthArray[DEPTHCOUNTER_ALL_R][displayIdx] += 1;
                  else
                    depthArray[DEPTHCOUNTER_ALL_F][displayIdx] += 1;

                  /* Find the match sequence base at this coord */
                  int sIdx = 0;
                  const char *seq = mspGetMatchSeq(msp);

                  if (mspGetMatchCoord(msp, qIdx, TRUE, numUnalignedBases, this, &sIdx))
                    {
                      /* Check we have the sequence. If not then don't do anything (this will
                       * show up as "unknown" in the read depth display) */
                      if (seq)
                        {
                          DepthCounter counter = getDepthCounterForChar(seq[sIdx - 1], msp->qStrand); // sIdx is 1-based

                          if (counter != DEPTHCOUNTER_NONE)
                            depthArray[counter][displayIdx] += 1;
                        }
                    }
                  else
                    {
                      /* No base here so it must be a gap in the match sequence */
                      if (msp->qStrand == BLXSTRAND_REVERSE)
                        depthArray[DEPTHCOUNTER_GAP_R][displayIdx] += 1;
                      else
                        depthArray[DEPTHCOUNTER_GAP_F][displayIdx] += 1;
                    }
                }
            }
        }
    } 
  
  /* Find the max and min depth (total depth over both strands) */
  minDepth = depthArray[DEPTHCOUNTER_ALL_F][0] + depthArray[DEPTHCOUNTER_ALL_R][0];
  maxDepth = minDepth;
  
  for (i = 1 ; i < displayLen; ++i)
    {
      const int cur_depth = depthArray[DEPTHCOUNTER_ALL_F][i] + depthArray[DEPTHCOUNTER_ALL_R][i];

      if (cur_depth < minDepth)
        minDepth = cur_depth;
      
      if (cur_depth > maxDepth)
        maxDepth = cur_depth;
    }  
}



/* Calculate the total depth of coverage of short-reads for the given range of ref seq coords.
 * depthArray must be the same length as displayRange. */
int BlxContext::calculateTotalDepth(const IntRange *range, const BlxStrand strand)
{
  int depth = 0;

  /* Loop through all MSP lists */
  for (int mspType = 0 ; mspType < BLXMSP_NUM_TYPES; ++mspType)
    {
      /* Only include MSPs of relevant types */
      if (!includeTypeInCoverage((BlxMspType)mspType))
        continue;
      
      /* Loop through all MSPs in this list */
      GArray *mspArray = featureLists[mspType];
    
      int i = 0;
      for (const MSP *msp = mspArrayIdx(mspArray, i); msp; msp = mspArrayIdx(mspArray, ++i))
        {
          /* If the alignment is in our range, increment the depth. Only include msps on the
           * given strand (or both strands if given strand is "none") */
          if ((strand == BLXSTRAND_NONE || msp->qStrand == strand) &&
              rangesOverlap(range, &msp->displayRange))
            ++depth;
        }
    } 

  return depth;
}

/* Utility to get the value from the depth array at the given coord for the given
 * counter. Validates the coord and counter are valid. The coord should be in display coords. */
int BlxContext::getDepthForCounter(const int coord, const DepthCounter counter)
{
  int result = 0;

  g_return_val_if_fail(valueWithinRange(coord, &fullDisplayRange) &&
                       counter > DEPTHCOUNTER_NONE &&
                       counter < DEPTHCOUNTER_NUM_ITEMS &&
                       depthArray[counter] != NULL,
                       result);

  int idx = invertCoord(coord, &fullDisplayRange, displayRev); // invert if display reversed
  idx -= fullDisplayRange.min(); // make 0-based

  result = depthArray[counter][idx];

  return result;
}


/* Return the read depth at the given display coord */
int BlxContext::getDepth(const int coord, 
                         const char *base_char,
                         const BlxStrand strand)
{
  int result = 0;
  g_return_val_if_fail(coord >= fullDisplayRange.min() &&
                       coord <= fullDisplayRange.max(), 
                       result);

  if (base_char && strand == BLXSTRAND_NONE)
    {
      /* Get the depth for the specific base for both strands */
      DepthCounter counter_f = getDepthCounterForChar(*base_char, BLXSTRAND_FORWARD);
      DepthCounter counter_r = getDepthCounterForChar(*base_char, BLXSTRAND_REVERSE);
      
      result = 
        getDepthForCounter(coord, counter_f) +
        getDepthForCounter(coord, counter_r);
    }
  else if (base_char && strand != BLXSTRAND_NONE)
    {
      /* Get the depth for the specific base and given strand */
      DepthCounter counter = getDepthCounterForChar(*base_char, strand);
      result = getDepthForCounter(coord, counter);
    }
  else if (strand == BLXSTRAND_NONE)
    {
      /* Get the depth for all reads for both strands */
      result = 
        getDepthForCounter(coord, DEPTHCOUNTER_ALL_F) + 
        getDepthForCounter(coord, DEPTHCOUNTER_ALL_R);
    }
  else if (strand == BLXSTRAND_FORWARD)
    {
      /* Get the depth for all reads for the forward strand */
      result = getDepthForCounter(coord, DEPTHCOUNTER_ALL_F);
    }
  else if (strand == BLXSTRAND_REVERSE)
    {
      /* Get the depth for all reads for the reverse strand */
      result = getDepthForCounter(coord, DEPTHCOUNTER_ALL_R);
    }
  else
    {
      /* All possible conditions should be covered above */
      g_warn_if_reached();
    }

  return result;
}


bool BlxContext::isSeqSelected(const BlxSequence *seq) const
{
  GList *foundItem = NULL;
  
  if (seq)
    {
      foundItem = g_list_find(selectedSeqs, seq);
    }
  
  return (foundItem != NULL);
}


/* Returns the group that the given sequence belongs to, if any (assumes the sequence
 * is only in one group; otherwise it just returns the first group it finds). */
SequenceGroup *BlxContext::getSequenceGroup(const BlxSequence *seqToFind) const
{
  SequenceGroup *result = NULL;
  
  if (!seqToFind)
    return result;
  
  /* Loop through all the groups until we find this sequence in one */
  GList *groupItem = sequenceGroups;
  for ( ; groupItem; groupItem = groupItem->next)
    {
      /* See if our sequence struct is in this group's list */
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      GList *foundItem = g_list_find(group->seqList, seqToFind);
      
      if (foundItem)
        {
          result = group;
          break;
        }
    }
  
  return result;
}


/* Return a list of all selected features of the given type. Result should be free'd by caller
 * using g_list_free */
GList *BlxContext::getSelectedSeqsByType(const BlxSequenceType type) const
{
  GList *result = NULL;

  GList *list_item = selectedSeqs;
  
  for ( ; list_item; list_item = list_item->next)
    {
      BlxSequence *curSeq = (BlxSequence*)(list_item->data);
      
      if (curSeq->type == type)
        {
          result = g_list_append(result, curSeq);
        }
    }

  return result;
}


/* If there is one (and only one) selected transcript then return it; otherwise return null. If
 * num_transcripts is given then return the number of selected transcripts. */
BlxSequence* BlxContext::getSelectedTranscript(int *num_transcripts) const
{
  BlxSequence *result = NULL;

  GList *list_item = selectedSeqs;
  int num_found = 0;
  
  for ( ; list_item; list_item = list_item->next)
    {
      BlxSequence *curSeq = (BlxSequence*)(list_item->data);
      
      if (curSeq->type == BLXSEQUENCE_TRANSCRIPT)
        {
          ++num_found;
          
          if (result)
            {
              /* Found more than one - don't know which to choose so return null */
              result = NULL;

              /* If we don't need to return the count, then exit now */
              if (!num_transcripts)
                break;
            }
          else
            {
              /* First one found: set the result. Continue to make sure there aren't any more */
              result = curSeq;
            }
        }
    }

  if (num_transcripts)
    *num_transcripts = num_found;

  return result;
}


void BlxContext::highlightBoxCalcBorders(GdkRectangle *drawingRect, 
                                         GdkRectangle *highlightRect,
                                         const IntRange *fullRange,
                                         const IntRange *highlightRange,
                                         const int yPadding)
{
  if (drawingRect && highlightRect && fullRange && highlightRange)
    {
      /* Get the full range in dna coords */
      IntRange fullDnaRange;
      convertDisplayRangeToDnaRange(fullRange, seqType, numFrames, displayRev, &refSeqRange, &fullDnaRange);
      
      /* Get the highlight range in dna coords */
      IntRange highlightDnaRange;
      convertDisplayRangeToDnaRange(highlightRange, seqType, numFrames, displayRev, &refSeqRange, &highlightDnaRange);
      
      /* Get the x coords for the start and end of the detail view display range */
      const int x1 = convertBaseIdxToRectPos(highlightDnaRange.min(true, displayRev), drawingRect, &fullDnaRange, TRUE, displayRev, TRUE);
      const int x2 = convertBaseIdxToRectPos(highlightDnaRange.max(true, displayRev), drawingRect, &fullDnaRange, TRUE, displayRev, TRUE);
      
      highlightRect->x = min(x1, x2);
      highlightRect->y = 0;

      highlightRect->width = abs(x1 - x2);

      highlightRect->height = 
        drawingRect->height + 
        roundNearest(charHeight() / 2.0) + 
        yPadding + 
        (2 * HIGHLIGHT_BOX_Y_PAD);
    }
}
