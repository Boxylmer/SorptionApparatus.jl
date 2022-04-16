
"""
    generatetemplate(apparatus; folder, filename)

Generate a blank apparatus Excel template.

* `folder` defaults to the working directory.
* `filename` defaults to a placeholder dependent on the apparatus type.
"""
function generatetemplate end

"""
    readtemplate(apparatus; filepath)

Read an apparatus template into its corresponding apparatus *setup*. 
"""
function readtemplate end

"""
    processtemplate(apparatus; template_filepath, results_filepath)
    processtemplate(apparatus; template_filepath)

Process an apparatus template directly from the template file. 
Returns a complete apparatus *system*, and if `results_filepath` is specified, will write the system to this path. 
"""
function processtemplate end