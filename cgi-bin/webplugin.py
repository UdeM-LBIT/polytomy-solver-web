#!/usr/bin/python
import os, sys, re, cgi, time

# settings.py
from settings import *

# Configure path accordingly
sys.path.insert(0, WEB_APP_FOLDER_PATH)

# python
from string import strip
from hashlib import md5
from pyvirtualdisplay import Display
import _mysql
import wsgiref.handlers
import json

# toolkits
from ete2 import WebTreeApplication, PhyloTree, faces, Tree, TreeStyle
from Bio.Phylo.Applications import _Phyml
from Bio.Phylo.TreeConstruction import _DistanceMatrix

# project
from algorithms.libpolytomysolver import PolytomySolver
from utils.TreeLib import TreeUtils, TreeClass

# start virtual display
xvfb=Display(visible=0, size=(1024, 768)).start()

# In order to extend the default WebTreeApplication, we define our own WSGI function that handles URL queries
def webplugin_app(environ, start_response, queries):
    asked_method = environ['PATH_INFO'].split("/")
    URL = environ['REQUEST_URI']

    # expected arguments from the URL (POST or GET method)
    treeid = queries.get("treeid", [None])[0]

    if not treeid:
        treeid = md5(str(time.time())).hexdigest()

    # NOTE : asked_methods are overriden by the base functions found in WebTreeApplication
    # (ex : catching "draw" here will not be executed)

    #
    # PolytomySolver
    #
    if asked_method[1]=="polytomysolver":
        speciesTree = queries.get("speciesTree", [None])[0]
        geneTree = queries.get("geneTree", [None])[0]
        geneSeq = queries.get("geneSeq", [None])[0]
        geneDistances = queries.get("geneDistances", [None])[0]
        sp_tol = queries.get("sp_tol", [None])[0]
        gn_ensembl_tree = queries.get("gn_ensembl", [None])[0]
        gn_reroot_mode = queries.get("gn_reroot_mode", [None])[0]
        gn_support_threshold = queries.get("gn_support_threshold", [None])[0]

        # Use the ensembl tree of life?
        if sp_tol != "0":
            f = open(WEB_APP_FOLDER_PATH+'/ressources/ensembl.nw', 'r')
            speciesTree = f.read()

        # Use an ensembl gene tree?
        if gn_ensembl_tree:
            geneTree = TreeUtils.fetch_ensembl_genetree_by_id(gn_ensembl_tree)

        # Contract low support branches (might want to add a y/n checkbox instead of running everytime with 0 as default)
        gn_tree_obj = TreeClass(TreeUtils.newick_preprocessing(geneTree)[0])
        gn_tree_obj.contract_tree(seuil=float(gn_support_threshold), feature='dist')

        # PolytomySolver v1.2.5
        # PolytomySolver(string speciesTreeString, string geneTreeString, string strDistances, string _rerootMode, bool _testEdgeRoots, bool _hasNonnegativeDistanceFlag, bool _useCache)
        if gn_reroot_mode in ["none","findbestroot","outputallroots"]:
            polytomysolver_out = PolytomySolver(str(speciesTree), gn_tree_obj.write(format=6).replace("%%",";;"), geneDistances, gn_reroot_mode, False, False, True)
        else:
            return '<b style="color:red;"> PolytomySolver reroot mode error. </b>'

        trees_processed = []

        # Get list of trees
        for line in polytomysolver_out.splitlines():
            if line[0] != "#":
                tree = TreeClass(line)
                speciesTree_obj = TreeClass(speciesTree)

                tree.set_species(sep="__",pos="postfix")
                tree.set_genes(sep="__",pos="prefix")

                lcamap = TreeUtils.lcaMapping(tree, speciesTree_obj)

                # Reconcile gene and species tree
                TreeUtils.reconcile(tree, lcamap)
                trees_processed.append(tree)
                if not application._load_tree(treeid, tree):
                    return '<b style="color:red;"> Cannot load the tree: %s' %treeid


        # Phyml - v20140520
        # Calculate the log likelihood of the output tree and the given gene sequence data
        # TODO : Wrap this into a separate function that will give log likelihood for a given tree and set of sequences
        if (geneSeq != None):
            # Write sequences to file
            with open("utils/phyml_tmp/%s.nex"%treeid, "w") as seqs_file:
                seqs_file.write(geneSeq)

            with open("utils/phyml_tmp/%s.newick"%treeid, "w") as newick_file:
                for tree in trees_processed:
                    # Write newick for trees in a .newick file named after the treeid of the first tree
                    # Convert all tree labels to gene names only
                    t_tmp = tree.copy()
                    leaves = t_tmp.get_leaves()
                    for leaf in leaves:
                        leaf.name=leaf.genes
                    newick_file.write(t_tmp.write(features=[])+"\n")

            # Set everything up to run phyml on the sequences and get the log likelihood for tree
            input_seqs = "%s/utils/phyml_tmp/%s.nex" %(WEB_APP_FOLDER_PATH, treeid)
            input_trees = "%s/utils/phyml_tmp/%s.newick" %(WEB_APP_FOLDER_PATH, treeid)

            phyml = _Phyml.PhymlCommandline(input=input_seqs, input_tree=input_trees, bootstrap=0)
	    phyml.set_parameter("-o","none")
            phyml.program_name = 'utils/phyml'
            phyml()

            # Fetch phyml output
            output_stats = "%s/utils/phyml_tmp/%s.nex_phyml_stats.txt" %(WEB_APP_FOLDER_PATH, treeid)
            output_tree = "%s/utils/phyml_tmp/%s.nex_phyml_tree.txt" %(WEB_APP_FOLDER_PATH, treeid)

            ll_keyword = ". Log-likelihood:"
            log_index = 0

            with open(output_stats) as search:
                for line in search:
                    if ". Log-likelihood:" in line:
                        line = line.replace(ll_keyword, "")
                        trees_processed[log_index].log_likelihood=line.strip()
                        trees_processed[log_index].tree_number=log_index+1
                        log_index += 1

            # Clean up tmp files
            os.remove(input_seqs)
            os.remove(input_trees)
            os.remove(output_stats)
            os.remove(output_tree)

            # Dropdown for all trees
            trees_dropdown ="""<script>
                            $("#select_trees_dropdown").on("change", function(){
                                var sel = document.getElementById("select_trees_dropdown");
                                var selected_value = sel.options[sel.selectedIndex].value;
                                var selected_value_json = JSON.parse(selected_value);
                                draw_tree(%s, selected_value_json.newick, "#img1", "name");
                                });

                                $(document).ready(function() {
                                $('#select_trees_dropdown').trigger("change");
                                });

                                </script>"""%treeid +\
                                        """<select id="select_trees_dropdown">"""

            trees_processed.sort(key=lambda x: x.log_likelihood)

            for tree in trees_processed:
                value = json.dumps({"newick":tree.write(features=[]),"log-likelihood":tree.log_likelihood})
                trees_dropdown += "<option value='%s'>Tree %d (Log-likelihood : %s)</option>"%(value,tree.tree_number,tree.log_likelihood)
            trees_dropdown += '</select>'

        return trees_dropdown


# ==============================================================================
# TREE LOADER
#
# This is my tree loading functions. I want the WebTreeApplication to
# use this method to load the trees
# ==============================================================================

# NOTE : Special consideration for current sample set... may not all be of this format!
def extract_species_code(name):
    name_lst = str(name).split("_")
    if len(name_lst) < 2:
        return name
    else:
        name_lst.pop(0)
        return ' '.join(name_lst).strip()

def my_tree_loader(tree):
    """ This function is used to load trees within the WebTreeApplication object. """

    t = PhyloTree(tree, sp_naming_function=None)

    #Check one leaf to see if species information is included
    if t.get_leaves()[0].species == "Unknown":
        t.set_species_naming_function(extract_species_code)

    return t

# ==============================================================================
# CUSTOM LAYOUTS
#
# These are my layout functions. I want the WebTreeApplication to use
# them for rendering trees
# ==============================================================================

LEAVE_FACES = [] # Global var that stores the faces that are rendered
                 # by the layout function
def main_layout(node):
    """ Main layout function. It controls what is shown in tree images. """

    # Add faces to leaf nodes. This allows me to add the faces from
    # the global variable LEAVE_FACES, which is set by the application
    # controler according to the arguments passed through the URL.
    if node.is_leaf():
        for f, fkey, pos in LEAVE_FACES:
            if hasattr(node, fkey):
                faces.add_face_to_node(f, node, column=pos, position="branch-right")
    else:
        # Add special faces on collapsed nodes
        if hasattr(node, "hide") and int(node.hide) == 1:
            node.img_style["draw_descendants"]= False
            collapsed_face = faces.TextFace(\
                    " %s collapsed leaves." %len(node), \
                    fsize=10, fgcolor="#444", ftype="Arial")
            faces.add_face_to_node(collapsed_face, node, 0)

    # Set node aspect. This controls which node features are used to
    # control the style of the tree. You can add or modify this
    # features, as well as their behaviour
    if node.is_leaf():
        node.img_style["shape"] = "square"
        node.img_style["size"] = 4
    else:
        node.img_style["size"] = 8
        node.img_style["shape"] = "sphere"

    # Evoltype: [1] Duplication A, [0] Speciation, [2] NAD or [-1] Losses.
    if hasattr(node,"type"):
        if str(node.type) == "2":
            node.img_style["fgcolor"] = "#FF8000"
            node.img_style["size"] = 12
        elif str(node.type) == "1":
            node.img_style["fgcolor"] = "#3F9933"
            node.img_style["size"] = 12
            node.img_style["hz_line_type"] = 0
            node.img_style["vt_line_type"] = 0
        elif str(node.type) == "0":
            node.img_style["fgcolor"] = "#000000"
            node.img_style["vt_line_color"] = "#000000"
            node.img_style["hz_line_color"] = "#000000"
        elif str(node.type) == "-1":
            node.img_style["fgcolor"] = "#FF0000"
            node.img_style["size"] = 12
            node.img_style["vt_line_color"] = "#FF0000"
            node.img_style["hz_line_color"] = "#FF0000"
            node.img_style["hz_line_type"] = 1
            node.img_style["vt_line_type"] = 1

    # If no evolutionary information, set a default style
    else:
        node.img_style["fgcolor"] = "#000000"
        node.img_style["vt_line_color"] = "#000000"
        node.img_style["hz_line_color"] = "#000000"
        node.img_style["hz_line_type"] = 0
        node.img_style["vt_line_type"] = 0

    # Parse node features features and conver them into styles. This
    # must be done like this, since current ete version does not allow
    # modifying style outside the layout function.
    if hasattr(node, "bsize"):
        node.img_style["size"]= int(node.bsize)

    if hasattr(node, "shape"):
        node.img_style["shape"]= node.shape

    if hasattr(node, "bgcolor"):
        node.img_style["bgcolor"]= node.bgcolor

    if hasattr(node, "fgcolor"):
        node.img_style["fgcolor"]= node.fgcolor


# ==============================================================================
# Checker function definitions:
#
# All checker actions must receive a node instance as unique argument
# and return True (node passes the filters) or False (node does not
# passes the filters).
#
# ==============================================================================

can_expand = lambda node: not node.is_leaf() and (hasattr(node, "hide") and node.hide==True)
can_collapse = lambda node: not node.is_leaf() and (not hasattr(node, "hide") or node.hide==False)
is_leaf = lambda node: node.is_leaf()
is_not_leaf = lambda node: not node.is_leaf()

# ==============================================================================
# Handler function definitions:
#
# All action handler functions must receive a node instance as unique
# argument. Returns are ignored.
#
# Note that there is a special action handler designed for searches
# within the tree. Handler receives node and searched term.
#
# ==============================================================================

def collapse(node):
    node.add_feature("hide", 1)
    node.add_feature("bsize", 25)
    node.add_feature("shape", "sphere")
    node.add_feature("fgcolor", "#bbbbbb")

def expand(node):
    try:
        node.del_feature("hide")
        node.del_feature("bsize")
        node.del_feature("shape")
        node.del_feature("fgcolor")
    except (KeyError, AttributeError):
        pass

def swap_branches(node):
    node.children.reverse()

def set_bg(node):
    node.add_feature("bgcolor", "#CEDBC4")


def set_as_root(node):
    node.get_tree_root().set_outgroup(node)

def phylomedb_clean_layout(node):
    phylomedb_layout(node)
    node.img_style["size"]=0

def search_by_feature(tree, search_term):
    ''' Special action '''
    attr, term = search_term.split("::")
    if not term:
        return None
    elif attr == "clean" and term == "clean":
        for n in tree.traverse():
            try:
                n.del_feature("bsize")
                n.del_feature("shape")
                n.del_feature("fgcolor")
            except:
                pass
    else:
        for n in tree.traverse():
            if hasattr(n, attr) and \
                    re.search(term,  str(getattr(n, attr)), re.IGNORECASE):
                        n.add_feature("bsize", 16)
                        n.add_feature("shape", "sphere")
                        n.add_feature("fgcolor", "#BB8C2B")


# ==============================================================================
# HTML generators
#
# Actions will be automatically added to the popup menus and attached
# to the action handler function. However, if just want to add
# informative items to the popup menu or external actions not
# associated to any handler, you can overwrite the default html
# generator of each action.
#
# html generators receive all information attached to the node and
# action in 5 arguments:
#
# * aindex: index of the action associated to this html generator
#
# * nodeid: id of the node to which action is attached
#
# * treeid: id of the tree in which node is present
#
# * text: the text string associated to the element that raised the
# action (only applicable to text faces actions)
#
# * node: node instance in which action will be executed.
#
#
# Html generator should return a text string encoding a html list
# item:
#
# Example: return "<li> my text </li>"
#
# ==============================================================================


def branch_info(aindex, nodeid, treeid, text, node):
    ''' It shows some info of the node in the popup menu '''
    return '''
           <li style="background:#eee; font-size:8pt;">
           <div style="text-align:left;font-weight:bold;">
            NODE ACTIONS
           </div>
            (<b>Branch: </b>%0.3f <b>Support:</b> %0.3f)<br>
           </li>'''  %\
                   (node.dist, node.support)

def search_in_ensmbl(aindex, nodeid, treeid, text, node):
    return '''<li>
              <a target="_blank" href="http://www.ensembl.org/common/Search/Results?species=all;idx=;q=%s">
              <img src=""> Search in ensembl: %s >
              </a>
              </li> ''' %\
                      (node.name, node.name)

def external_links_divider(aindex, nodeid, treeid, text, node):
    ''' Used to show a separator in the popup menu'''
    if node.is_leaf():
        return '''<li
        style="background:#eee;font-size:8pt;"><b>External
        links</b></li>'''
    else:
        return ""

def topology_action_divider(aindex, nodeid, treeid, text, node):
    return '''<li style="background:#eee;"><b>Tree node actions</b></li>'''

# ==============================================================================
# TREE RENDERER
#
# By default, ETE will render the tree as a png image and will return
# a simplistic HTML code to show the image and interact with
# it. However, it is possible to wrap such functionality to preprocess
# trees in a particular way, read extra parameters from the URL query
# and/or produce enriched HTML applications.
#
# Tree renderer wrappers receive the tree object, its id, and the WSGI
# application object. They MUST call the
# application._get_tree_img(tree) method and return the a HTML
# string.
#
# A simplistic wrapper that emulates the default WebTreeApplication
# behaviour would be:
#
# def tree_renderer(tree, treeid, application):
#    html = application._get_tree_img(treeid = treeid)
#    return html
#
# ==============================================================================
def tree_renderer(tree, treeid, application):
    # The following part controls the features that are attached to
    # leaf nodes and that will be shown in the tree image. Node styles
    # are set here, and faces are also created. The idea is the
    # following: user can pass feature names using the URL argument
    # "tree_features". If the feature is handled by our function and
    # it is available in nodes, a face will be created and added to
    # the global variable LEAVE_FACES. Remember that our layout
    # function uses such variable to add faces to nodes during
    # rendering.

    # Extracts from URL query the features that must be drawn in the tree
    asked_features = application.queries.get("show_features", ["name"])[0].split(",")
    print >> sys.stderr, asked_features

    def update_features_avail(feature_key, name, col, fsize, fcolor, prefix, suffix):
        text_features_avail.setdefault(feature_key, [name, 0, col, fsize, fcolor, prefix, suffix])
        text_features_avail[feature_key][1] += 1

    tree.add_feature("fgcolor", "#000000")
    tree.add_feature("shape", "sphere")
    tree.add_feature("bsize", "8")
    tree.dist = 0

    # This are the features that I wanto to convert into image
    # faces. I use an automatic function to do this. Each element in
    # the dictionary is a list that contains the information about how
    # to create a textFace with the feature.
    leaves = tree.get_leaves()
    formated_features =  {
            # feature_name: ["Description", face column position, text_size, color, text_prefix, text_suffix ]
            "name": ["name", len(leaves), 0, 11, "#000000", "", ""],
            }

    # populates the global LEAVE_FACES variable
    global LEAVE_FACES
    LEAVE_FACES = []
    unknown_faces_pos = 2
    for fkey in asked_features:
        if fkey in formated_features:
            name, total, pos, size, color, prefix, suffix = formated_features[fkey]
            f = faces.AttrFace(fkey, ftype="Arial", fsize=size, fgcolor=color, text_prefix=prefix, text_suffix=suffix)
            LEAVE_FACES.append([f, fkey, pos])
        else:
            # If the feature has no associated format, let's add something standard
            prefix = " %s: " %fkey
            suffix = ","
            if fkey == asked_features[-1]:
                suffix=""
            f = faces.AttrFace(fkey, ftype="Arial", fsize=10, fgcolor="#666666", text_prefix=prefix, text_suffix=suffix)
            LEAVE_FACES.append([f, fkey, unknown_faces_pos])
            unknown_faces_pos += 1

    text_features_avail = {}
    for l in leaves:
        for f in l.features:
            if not f.startswith("_"):
                text_features_avail.setdefault(f, 0)
                text_features_avail[f] = text_features_avail[f] + 1

    html_features = '''
      <div id="tree_features_box">
      <div class="tree_box_header">Available tree features
      <img src="webplugin/close.png" onclick='$(this).closest("#tree_features_box").hide();'>
      </div>
      <form id="form_tree_features" action='javascript: set_tree_features("", "", "");'>

      '''

    for fkey, counter in text_features_avail.iteritems():
        if fkey in asked_features:
            tag = "CHECKED"
        else:
            tag = ""

        fname = formated_features.get(fkey, [fkey])[0]

        html_features += '<INPUT NAME="tree_feature_selector" TYPE=CHECKBOX %s VALUE="%s">%s (%s/%s leaves)</input><br> ' %\
                (tag, fkey, fname, counter, len(leaves))


    html_features += '''<input type="submit" value="Refresh"
                        onclick='javascript:
                                // This piece of js code extracts the checked features from menu and redraw the tree sending such information
                                var allVals = [];
                                $(this).parent().children("input[name=tree_feature_selector]").each(function(){
                                if ($(this).is(":checked")){
                                    allVals.push($(this).val());
                                }});
                                draw_tree("%s", "", "#img1", allVals.join(",") );'
                       >
                       </form></div>''' %treeid

    features_button = '''
     <li><a href="javascript:;" onclick='show_box(event, $(this).closest("#tree_panel").children("#tree_features_box"));'>
     <img width=16 height=16 src="webplugin/icon_tools.png" alt="Select Tree features">
     </a></li>'''

    download_button = '''
     <li><a href="webplugin/tmp/%s.png" target="_blank">
     <img width=16 height=16 src="webplugin/icon_attachment.png" alt="Download tree image">
     </a></li>''' %treeid

    search_button = '''
      <li><a href="javascript:;" onclick='javascript:
          var box = $(this).closest("#tree_panel").children("#search_in_tree_box");
          show_box(event, box); '>
      <img width=16 height=16 src="webplugin/icon_search.png" alt="Search in tree">
      </a></li>'''

    clean_search_button = '''
      <li><a href="javascript:;" onclick='run_action("%s", "", %s, "clean::clean");'>
      <img width=16 height=16 src="webplugin/icon_cancel_search.png" alt="Clear search results">
      </a></li>''' %\
              (treeid, 0)


    search_select = '<select id="ete_search_target">'
    # Name first
    if 'name' in text_features_avail:
        search_select += '<option value="name">name</option>'
        del text_features_avail['name']
    for fkey in text_features_avail:
        search_select += '<option value="%s">%s</option>' %(fkey,fkey)
    search_select += '</select>'

    main_search = '''
    <form onsubmit='javascript:
                     search_in_tree("%s", "%s",
                                    $(this).closest("form").children("#ete_search_term").val(),
                                    $(this).closest("form").children("#ete_search_target").val());'
          action="javascript:void(0);">
          <input id="ete_search_term" type="text" value="Search"
          style="font-style:italic;color:grey;"
          onfocus="if(this.value == 'Search') { this.value = ''; this.style.color='black'; this.style.fontStyle='normal'}"'>%s</input>
     </form>
     '''%\
             (treeid, 0, search_select)


    search_form = '''
     <div id="search_in_tree_box">
     <div class="tree_box_header"> Search in Tree
     <img src="webplugin/close.png" onclick='$(this).closest("#search_in_tree_box").hide();'>
     </div>
     <form onsubmit='javascript:
                     search_in_tree("%s", "%s",
                                    $(this).closest("form").children("#ete_search_term").val(),
                                    $(this).closest("form").children("#ete_search_target").val());'
          action="javascript:void(0);">
     %s
     <input id="ete_search_term" type="text" value=""'></input>
     <br><i>(Searches are not case sensitive and accept Perl regular expressions)</i>
     <br>
     </form>
     <i> (Press ENTER to initiate the search)</i>
     </div>
     ''' %\
             (treeid, 0, search_select) # 0 is the action index associated to the
                                        # search functionality. This means that this action is the
                                        # first to be registered in WebApplication.


    buttons = '<div id="ete_tree_buttons">' +\
            main_search + features_button + clean_search_button + download_button +\
            '</div>'

    tree_panel_html = '<div id="tree_panel">' + search_form + html_features + buttons + '</div>'

    # Now we render the tree into image and get the HTML that handles it
    tree_html = application._get_tree_img(treeid = treeid)

    newick = '''<br><textarea  id="newick_out" title="Newick" style="width: 500px; height: 80px; border:dashed 1px;">''' + \
            tree.write(features=[]) +\
            "</textarea>"

    legend = '''<br>
    <ul>
    <il><img src="webplugin/legend/spec.png"/>Speciation</li>
    <il><img src="webplugin/legend/dup.png"/>Duplication</li>
    <il><img src="webplugin/legend/loss.png"/>Loss</li>
    <il><img src="webplugin/legend/nad.png"/>Non-Apparent Duplication</li>
    </ul>
    '''

    try:
        tree_ll = '''<br> Log-likelihood : '''+ tree.log_likelihood
    except Exception as e:
        tree_ll = ""

    # Now that the tree image (and final html) have been generated, stop the virtual X server.
    xvfb.stop()

    # Let's return enriched HTML
    return tree_panel_html + tree_html + legend + newick + tree_ll

# ==============================================================================
#
# Main WSGI Application
#
# ==============================================================================
if __name__ == '__main__':

	# Create a basic ETE WebTreeApplication
	application = WebTreeApplication()

	# Set your temporal dir to allow web user to generate files. This two
	# paths should point to the same place, one using the absolute path in
	# your system, and the other the URL to access the same
	# directory. Note that the referred directory must be writable by the
	# webserver.
	application.CONFIG["temp_dir"] = TMP_DIR_PATH
	application.CONFIG["temp_url"] = TMP_WEB_RELATIVE_PATH # Relative to web site Document Root

	# Set the DISPLAY port that ETE should use to draw pictures. You will
	# need a X server installed in your server and allow webserver user to
	# access the display. If the X server is started by a different user
	# and www-data (usally the apache user) cannot access display, try
	# modifiying DISPLAY permisions by executing "xhost +"
	# Alternatively, use pyvirtualdisplay for headless X (as used here).
	application.CONFIG["DISPLAY"] = str(xvfb.new_display_var)

	# We extend the minimum WebTreeApplication with our own WSGI
	# application
	application.set_external_app_handler(webplugin_app)

	# Lets now apply our custom tree loader function to the main
	# application
	application.set_tree_loader(my_tree_loader)

        # Set tree style
        ts = TreeStyle()
        ts.show_scale = False
	application.set_tree_style(ts)

	# And our layout as the default one to render trees
	application.set_default_layout_fn(main_layout)

	# I want to make up how tree image in shown using a custrom tree
	# renderer that adds much more HTML code
	application.set_external_tree_renderer(tree_renderer)

	# ==============================================================================
	# ADD CUSTOM ACTIONS TO THE APPLICATION
	#
	# The function "register_action" allows to attach functionality to
	# nodes in the image. All registered accions will be shown in the
	# popup menu bound to the nodes and faces in the web image.
	#
	#
	# register_action(action_name, target_type=["node"|"face"|"layout"|"search"], \
	#                 action_handler, action_checker, html_generator_handler)
	#
	# When the Application is executed it will read your registered
	# acctions and will do the following:
	#
	# 1. Load the tree and get the image map
	#
	# 2. For each node and face in the tree, it will browse all registered
	# actions and will run the action_checker function to determine if the
	# action must be activated for such node or face
	#
	# 3. If action_checker(node) returns True, the action will be attached
	# to the context menu of that specific node or face, otherwise it will
	# be hidden.
	#
	# 4. When a click is done on a specific node, popup menus will be
	# built using their active actions. For this, ETE will use the
	# html_generator function associated to each function if
	# available. Otherwise, a popup entry will be added automatically.
	#
	# 5. When a certain action is pressed in the popup menus, the
	# action_handler function attached to the action will be executed over
	# its corresponding node, and the tree image will be refreshed.
	#
	# Special values:
	#
	#  action_checker = None : It will be interpreted as "Show allways"
	#  html_generator = None : Autogenerate html and link to action
	#  action_handler = None : Action will be ignored
	#
	# ==============================================================================

	# We first register the special action "search" which is attached to
	# our custom search function.
	application.register_action("", "search", search_by_feature, None, None)

	# Node manipulation options (bound to node items and all their faces)
	application.register_action("branch_info", "node", None, None, branch_info)
	application.register_action("Highlight background", "node", set_bg, None, None)
	application.register_action("Swap children", "node", swap_branches, is_not_leaf, None)
	#application.register_action("Set as root", "node", set_as_root, None, None)
	#application.register_action("<b>Collapse</b>", "node", collapse, can_collapse, None)
	#application.register_action("Expand", "node", expand, can_expand, None)

	# Actions attached to node's content (shown as text faces)
	#application.register_action("divider", "face", None, None, external_links_divider)

	#application.register_action("Default layout", "layout", main_layout, None, None)
	#application.register_action("Clean layout", "layout", main_layout, None, None)

	wsgiref.handlers.CGIHandler().run(application)
