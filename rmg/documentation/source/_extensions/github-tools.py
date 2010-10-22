"""
:Description: Sphinx extension to remove leading under-scores from directories names in the html build output directory.
"""
import os
import shutil


def setup(app):
    """
    Add a html-page-context  and a build-finished event handlers
    """
    app.connect('html-page-context', change_pathto)
    app.connect('build-finished', move_private_folders)
                
def change_pathto(app, pagename, templatename, context, doctree):
    """
    Replace pathto helper to change paths to folders with a leading underscore.
    """
    pathto = context.get('pathto')
    def gh_pathto(otheruri, *args, **kw):
        if otheruri.startswith('_'):
            otheruri = otheruri[1:]
        return pathto(otheruri, *args, **kw)
    context['pathto'] = gh_pathto
    if context.has_key('body'):
        context['body'] = context['body'].replace(u'_images/',u'images/')
    
def move_private_folders(app, e):
    """
    remove leading underscore from folders in in the output folder.
    
    Only does anything if output folder is 'html'
    """
    import shutil
    junk, outdirname = os.path.split(app.builder.outdir)
    if outdirname != 'html': return
    def join(dir):
        return os.path.join(app.builder.outdir, dir)
    for item in os.listdir(app.builder.outdir):
        if item.startswith('_') and os.path.isdir(join(item)):
            if os.path.exists(join(item[1:])):
                shutil.rmtree(join(item[1:]))
            shutil.move(join(item), join(item[1:]))
