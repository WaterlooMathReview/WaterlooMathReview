CIS-Dept. Lean Skin for PmWiki 2
================================
README.txt

About the skin
--------------
This is a clean, minimalist skin for PmWiki 2 that was created with usability
in mind. It evolved from a similar skin for PmWiki 1, which can be found here:
   http://www.pmwiki.org/wiki-v1/Cookbook/LeanSkin

Design objectives
-----------------
* Be understated without being stripped bare.  Emphasize the content.
* Treat screen pixels as valuable and scrolling as inconvenient.
* Avoid displaying nonessential duplicate links.
* Avoid linking from a page to itself (self-referencing links).
* Give navigation links conspicuous placement in the layout.
* Be easy to use for both inexperienced and seasoned veteran wiki users.
* Allow wiki links and last-modified output to be obscure or invisible.
* Provide an adequate text area (approx. 80 x 30) for easy editing.
* Keep the main content area a reasonable width for easier reading.
* Produce pages that render similarly in various browsers.
* Encourage wiki users to create well-structured web documents.
* Be easily customizable to suit personal preference.

The skin generates valid XHTML 1.0 and CSS output according to w3c.org's
validation services at http://validator.w3.org/ and
http://jigsaw.w3.org/css-validator/ respectively.

Access Keys
-----------
The skin includes these access keys:

accesskey+'c' = Recent [C]hanges in current group 
accesskey+'a' = [A]ll recent changes
accesskey+'b' = Edit side[B]ar
accesskey+'e' = [E]dit Page
accesskey+'h' = Page [H]istory
accesskey+'p' = [P]review page when editing
accesskey+'s' = [S]ave page when editing

Using access keys allows your site to have "Invisible Wiki Links",
where the wiki links aren't displayed but their access keys still work
(in Firefox but not Internet Explorer). Just uncomment the appropriate
lines in the stylesheet (lean.css) and you will have a quick-and-
easy-to-edit wiki that looks like an ordinary web site.


Installing the skin
-------------------
See INSTALL.txt


Revision History
----------------
Ver. 0.8
* The initial release of the skin for PmWiki 2.

Ver. 0.9
* Stylesheet settings can now override PmWiki defaults.
* Switched to $SkinDirUrl instead of $SkinDir.
* Now SideBar bold text has a margin for better appearance.
 
Ver. 0.9.1
* Changed site-related links to non-bold text.

Ver. 0.10
* Added "media print {}" block for printing from browse view.
* Improved wrapping (in Firefox) of History Page diff code.
* Added a copyright (/copyleft) notice capability.
* Switched to $DefaultPage from Main.HomePage to in lean.php.
* Integrated "Credit Block" PmWiki add-on.
* Fixed bug where $SkinTextHomeLink </div> was missing its '/'.
* Now supporting some <!--Page...Fmt--> directives, so
  (:noheader:) and (:nofooter:) special markups will work.
* Removed print/print.php script, since it wasn't being run. 

Ver 0.11
* Improved and refined the stylesheet, with lots of testing.
* Added a capability to use smaller fonts.
* Added a capability to use non-bold, underlined hyperlinks.
* The new default font size may seem tiny compared to some of the
  other skins in the PmWiki Cookbook, but not compared to the fonts
  of high-traffic web sites I visited (Amazon, eBay, Google, CNN,
  Yahoo, Wikipedia, and a bunch of others) in order to determine
  an ideal default size. 

Ver. 0.12
* Switched the font in the main content area to Verdana. Now the text
  (particularly small text) is easier to read and the same screen area
  holds more content. 
* The editable text area of the Edit page is shorter by 20px for less
  scrolling on 1024x768 displays. Chose a mid-sized font size as default
  for the main content page. Smaller and larger optional sizes are available.
* Fine-tuned headings' sizes for each of the content-area font sizes.
* Now providing a link back to a page when editing or viewing history.

Ver. 0.12.1
* Reworked the Printable View skin.

Ver. 0.13
* Now suppressing the header and SideBar on the Edit and History pages.
  Expanded the size of the Edit page's text area.
* Made some changes to encourage wiki pages to be well-structured web
  documents (where there is one <h1> heading per page and the heading
  levels are used to indicate sections' importance in the document).
  - The page title in the header area is now a <h1> heading.
  - SideBar headings can now be <h6> or <h1> headings.
* Vertical space above and below headings is now specified in the
  stylesheet to (enhance document appearance and improve cross-browser
  consistency).
* Added unique identifiers (id='...') for some items in the header and
  footer.  Now those items' appearance may be controlled (e.g. displayed
  or not) via the CSS stylesheet.
* Browser Compatibility:
  - Works in Firefox 1.0 / Win98 SE.
  - Works in IE 5.5 / Win98 SE.
  - Works in Firefox 1.0 / Debian Woody+Backports.org.
  - Works in Konqueror 2.2.2 / Debian Woody+Backports.org.

Ver. 0.14
* The overall width is slightly wider.  The extra width was added to the
  wikitext area.
* Diff output on History Pages is wider (same width as textarea on Edit
  Pages now).
* Preview output on Edit Pages is narrower (same width as wikitext area now).
* Heading margins are smaller.
* History Page links to Edit Page and vice versa.

Ver. 0.14.1
* Adjusted CSS so wiki style markup can more readily override the stylesheet.

Ver. 0.14.2
* Fixed IE compatibility bug introduced in v0.14.1.

Ver. 0.14.3
*  Improved links on Edit and History Pages.  (Suggested by Jeff Schuler.)

Ver. 0.15
* Adapted the skin so it works on fresh installations of PmWiki 2.0.beta44
  or newer, which use the Site group rather than the Main group for site-
  related pages.

Ver. 0.16
* Added a searchbox to the page header.
* Added CSS for compatibility with v2.0.0 SideBar wikistyle.

Ver. 0.16.1
* Added preview button to the edit form.
* Adjusted fonts in the stylesheet.

Ver. 0.16.2
* Shortcut icon (favicon) is now configurable via $EnableSkinIcon setting.
* Trail path separator is ' &gt; ' by default (instead of PmWiki's ' | ').
* $Action has been replaced with $ActionTitle to restore correct title.
* CSS stylesheet rearranged so it's now overridable by the administrator.

Ver. 0.16.3
* Using a logo or text for the home link is configurable with $EnableSkinLogo.
* The URL of the logo image file is configurable with $LeanLogoFile.


Author
------
Hagan Fox - http://Qdig.SourceForge.net  (The site will look familiar. ;) )

