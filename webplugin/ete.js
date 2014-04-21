/*  requires jquery */

var ete_webplugin_URL = "http://127.0.0.1:8080/wsgi/webplugin.py";
var loading_img = '<img border=0 src="/webplugin/loader.gif">';

function draw_tree(treeid, newick, recipient, extra_params){
    var params = {"tree": newick, "treeid": treeid};
    if ( extra_params != undefined ){
        var params =  $.extend(params, extra_params);
    }
    $(recipient).html(loading_img);
    $(recipient).load(ete_webplugin_URL+'/draw', params);
}

function polytomysolver(treeid, speciesTree, geneTree, distances, show_features, sp_tol, sp_ensembl, recipient, extra_params){
    var params = {"speciesTree": speciesTree,
        "geneTree": geneTree,
        "distances": distances,
        "treeid": treeid,
        "show_features": show_features,
        "sp_tol":sp_tol,
        "sp_ensembl": sp_ensembl};

    if ( extra_params != undefined ){
        var params =  $.extend(params, extra_params);
    }

    $(recipient).fadeIn().html(loading_img);

    $(recipient).load(ete_webplugin_URL+'/polytomysolver', params,

            function( response, status, xhr ) {

                if ( status == "error" ) {
                    $( recipient ).html( '<b style="color:red;"> Oops! Something went wrong.</b>' );
                }else{
                    $(recipient).css("display","none");
                    $(recipient).find("img.ete_tree_img").load(function(){
                        $(recipient).fadeIn('short').css({"display":"inline-block" ,"min-height":$(recipient).height(),"height":"auto"});
                    });
                }});
}

function show_context_menu(treeid, nodeid, actions, textface){
    if ( textface==undefined ){
        var textface = "";
    }
    var params = {"treeid": treeid, "show_actions": actions, "nid": nodeid, "textface": textface};
    $("#popup").html(loading_img);
    $('#popup').load(ete_webplugin_URL+'/get_menu', params);
}

function run_action(treeid, nodeid, aindex, search_term){
    var recipient = "#ETE_tree_"+treeid;

    $(recipient).css({"min-height":$(recipient).height()});
    $(recipient).fadeIn().html(loading_img);

    var params = {"treeid": treeid, "nid": nodeid, "aindex": aindex, "search_term": search_term};
    $(recipient).fadeOut().fadeIn('short').load(ete_webplugin_URL+'/action', params,
            function( response, status, xhr ) {
                if ( status == "error" ) {
                    $( recipient ).html( '<b style="color:red;"> Oops! Something went wrong.</b>' );
                } else {
                    // $(recipient).fadeIn('short');
                }
            });
}

function random_tid(){
    return Math.ceil(Math.random()*10000000);
}

function bind_popup(){
    $(".ete_tree_img").bind('click',function(e){
        $("#popup").css('left',e.pageX-2 );
        $("#popup").css('top',e.pageY-2 );
        $("#popup").css('position',"absolute" );
        $("#popup").css('background-color',"#fff" );
        $("#popup").draggable({ cancel: 'span,li' });
        $("#popup").show();
    });
}

function hide_popup(){
    $('#popup').hide();
}

function search_in_tree(treeid, search_index_action, search_term, term_target){
    var term = term_target + "::" + search_term;
    run_action(treeid, "", search_index_action, term);
}

function show_box(e, box) {
    box.draggable({ handle: '.text_header_box' });
    box.css('left',e.pageX+5 );
    box.css('top',e.pageY+5 );
    box.css('position',"absolute" );
    box.show();
}

$(document).ready(function(){
    hide_popup();
});