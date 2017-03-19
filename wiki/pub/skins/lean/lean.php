<?php if (!defined('PmWiki')) exit();
/**
* This is lean.php, the php script portion of the Lean Skin for PmWiki 2.
*/

## Add a custom page storage location for the
## custom Edit Form and a Preferences page.
global $WikiLibDirs;
$PageStorePath = dirname(__FILE__)."/wikilib.d/\$FullName";
$where = count($WikiLibDirs);
if ($where>1) $where--;
array_splice($WikiLibDirs, $where, 0,
  array(new PageStore($PageStorePath)));

## Enable the Preferences page.
global $PageEditForm, $XLLangs;
XLPage('light', 'Site.LeanXLPage');
array_splice($XLLangs, -1, 0, array_shift($XLLangs));

## Enable the skin's custom EditForm, either
## configurable via the prefs page or not.
$EnableEditFormPrefs = TRUE;
global $PageEditForm, $XLLangs;
if ($EnableEditFormPrefs == TRUE) {
  SDV($PageEditForm, '$[Site.EditForm]');
} else {
  SDV($PageEditForm, 'Site.LeanEditForm');
}

## For backward compatibility.
global $ActionTitle, $Action; if (!$ActionTitle) { $ActionTitle = $Action; }

## Favicon
global $EnableSkinIcon;
SDV($EnableSkinIcon, 1);
if ($EnableSkinIcon == TRUE) {
  global $SkinIconFmt; $SkinIconFmt = "
  <link rel='icon' href='\$SkinDirUrl/leanicon.gif' type='image/gif' />
  <link rel='SHORTCUT ICON' href='\$SkinDirUrl/leanicon.gif' />"; }
    
## Trail path separator
global $TrailPathSep;  SDV($TrailPathSep, ' &gt; ');

global $HTMLHeaderFmt; $HTMLHeaderFmt['leanskin'] =
   "  <link rel='stylesheet' href='\$SkinDirUrl/lean.css' type='text/css' />\n  ";

## Use GUI buttons on edit pages, including add some extra buttons.
global $EnableGUIButtons, $GUIButtons;
$EnableGUIButtons = 1;
$GUIButtons['h3'] = array(402, '\\n!!! ', '\\n', '$[Subheading]',
                     '$GUIButtonDirUrlFmt/h3.gif"$[Subheading]"');
$GUIButtons['indent'] = array(500, '\\n->', '\\n', '$[Indented text]',
                     '$GUIButtonDirUrlFmt/indent.gif"$[Indented text]"');
$GUIButtons['ul'] = array(530, '\\n* ', '\\n', '$[Unordered list]',
                     '$GUIButtonDirUrlFmt/ul.gif"$[Unordered (bullet) list]"');
$GUIButtons['ol'] = array(520, '\\n# ', '\\n', '$[Ordered list]',
                     '$GUIButtonDirUrlFmt/ol.gif"$[Ordered (numbered) list]"');
$GUIButtons['table'] = array(600,
                     '(:table border=1 width=80%:)\\n(:cell style=\'padding:5px\;\':)\\n1a\\n(:cell style=\'padding:5px\;\':)\\n1b\\n(:cellnr style=\'padding:5px\;\':)\\n2a\\n(:cell style=\'padding:5px\;\':)\\n2b\\n(:tableend:)\\n', '', '',
                     '$GUIButtonDirUrlFmt/table.gif"$[Table]"');
#$GUIButtons['table'] = array(600,
#                     '||border=1 width=80%\\n||!Hdr ||!Hdr ||!Hdr ||\\n||     ||     ||     ||\\n||     ||     ||     ||\\n', '', '',
#                     '$GUIButtonDirUrlFmt/table.gif"$[Table]"');

## Add a (:noleft:) directive.
Markup('noleft','directives','/\\(:noleft:\\)/e',
  "PZZ(\$GLOBALS['PageLeftFmt']='')");

## Add a (:noright:) directive.
Markup('noright','directives','/\\(:noright:\\)/e',
  "PZZ(\$GLOBALS['PageRightFmt']='')");

# Link back to the page from Edit and History pages.
global $SkinPageLinkPreFmt, $SkinPageLinkPostFmt, $SkinHideHead, $SkinHideSide,
  $SkinPageTitlePre, $SkinPageTitlePost;
if (@in_array($_GET['action'], array('edit', 'diff'))) {
    $SkinPageTitlePre    = '';
    $SkinPageTitlePost   = '';
    $SkinPageLinkPreFmt  = "$[Return to] <a href='\$PageUrl'>";
    $SkinPageLinkPostFmt = "</a> &nbsp;(<a
       href='\$PageUrl?action=edit'>$[Edit]</a>)";
    $SkinHideHead        = "style='display:none;'";
    $SkinHideSide        = "style='display:none;'";
} else {
    $SkinPageTitlePre    = '<h1>';
    $SkinPageTitlePost   = '</h1>';
    $SkinPageLinkPreFmt  = "<span style='display:none;'>";
    $SkinPageLinkPostFmt = '</span>';
    $SkinHideHead        = '';
    $SkinHideSide        = '';
}
if (@$_GET['action'] == 'edit' ) {
    $SkinPageLinkPostFmt  = "</a> &nbsp;(<a
       href='\$PageUrl?action=diff'>$[History]</a>)";
}

# Don't link to the Home Page when we're on it.
# Works for both a logo link ($SkinLogoHomeLink) or a text link ($SkinTextHomeLink).
global $DefaultPage, $ScriptUrl, $SkinHomeLink, $WikiTitle, $EnableSkinLogo, $LeanLogoFile;
SDV($EnableSkinLogo, 0);
SDV($LeanLogoFile, "$SkinDirUrl/leanlogo.gif");
if (empty($pagename) || $pagename==$DefaultPage) {
  $SkinLogoHomeLink = "<div id='sitelogo'><img src='\$LeanLogoFile'
        alt='$WikiTitle' border='0' /></div>";
  $SkinTextHomeLink = "<div id='sitelogotext'>$WikiTitle</div>";
} else {
  $SkinLogoHomeLink = "<div id='sitelogo'><a href='\$ScriptUrl'><img
        src='\$LeanLogoFile'
        alt='$WikiTitle' border='0' /></a><div id='sitelogo'></div>";
  $SkinTextHomeLink = "<div id='sitelogotext'><a href='\$ScriptUrl'>$WikiTitle</a></div>";
}
if ($EnableSkinLogo == TRUE) {
  $SkinHomeLink = $SkinLogoHomeLink;
} else {
  $SkinHomeLink = $SkinTextHomeLink;
}

## Don't link to the Group.Group page if it's already the curren page.
global $ScriptUrl, $GroupHdrFmt;
$page_array = explode('.',$pagename);
if (empty($pagename)
  || $pagename == $DefaultPage
  || $page_array['0'] == $page_array['1']
  || $page_array['1'] == 'HomePage') {
 $GroupHdrFmt = "\$Groupspaced";
} else {
 $GroupHdrFmt = "<a href='\$ScriptUrl/\$Group'
          title='\$Groupspaced \$[Home]'>\$Groupspaced</a>";
}

# Copyright notice
# TODO Language $[settings] don't work
global $SkinCopyright;
/*
$Copyright = "All text is available under the terms of the
          <a href='http://www.gnu.org/copyleft/fdl.html'
           title='GNU FDL Home'>GNU Free Documentation License</a>";
*/

if (isset($Copyright)) {
  $SkinCopyright = "<span id='copyright' title='Copyright notice'>$Copyright</span>";
} else {
  $SkinCopyright = '';
}

?>
