# BASH scripting

Bash scripts are a nice, basic way to automate, record your commands, and start practicing reproducibility!

## Basics

A bash script is really just a series of bash commands and/or any other command that's been configured to run via the shell. They are executed in the order you put them in the script. The very first line should be the "shebang" `#!`, telling the computer which interpreter to use (this is only really necessary if you want to use a different interpreter/version/whatever than the default, but put it in there anyway).

Shebang:

    #!/bin/bash
    
You need a text editor (discussed below) to write the script. For good practice and human readability, the file should be saved with the `.sh` extension.

After the shebang, any line that begins with a `#` is considered a comment and is ignored by the computer when you run the script. 

Scripts are run as follows:

    $ bash script.sh
    
### Text editors

Anything that will save your file as a plain text. There are several options on the server already (nano, emacs, vim). These are applications to be used in the terminal *only* and the files are saved in the current directory you're in when you open the application. On the other hand, there are several options orf GUIs (TextWrangler, XCode, TextEdit, etc) that you can use on the local computer. These feel more intuitive to people, but you must then copy the script onto the server (using `scp`) if you want to execute it there. 

### Variables

    varName=value

- varName stores "value" for the duration of your shell session
- No spaces around =
- varName can have letters, numbers and underscores but cannot start with a number
- Retrieve the value by prepending variable name with $ (ie $varName)
- $ can be used with curly braces or double quotes to avoid confusion with any following text (ie ${varName})
- These are super useful at the beginning of a script to store values for arguments that you might want to change when you run the script on different data.

Examples (run in order):

    $ bio373=/sratch/bio373_2020
    $ workdir=$bio373/data
    $ echo $workdir
    $ ls ${workdir}
    $ test=AAA
    $ test_var=BBB
    $ echo $test_var
    $ echo ${test}_var

### Control the flow of your commands

You can apply basic programming principles in bash scripts. Bash has it's own syntax and the following are some examples of the bash syntax.

---------------------------------------

**Conditional**

A space between the square bracket and the conditional statement is required!

    if [ <condition> ]
    then
        command
    else
        otherCommand
    fi

Example:

    if [ ! -d $outdir ]
    then
        mkdir ${outdir}
    fi
    
If the path stored in the $outdir variable does not exist, make it.

---------------------------------------

**Loops**

    for item in list
    do
        command $item
    done
    
Example:

    for genotype in "w51690" "w66039"
    do
        echo "Genotype $genotype has these files in this directory: "
        ls ${genotype}*
    done
    
Another one I like to use takes a file as input and loops through each line of the file.

    while read item
    do
        command ${item}
    done < input.txt    

Example:

    while read accession
	    do
	        echo ${accession} 
	        ~/bin/fastq-dump --origfmt -I --split-files --gzip -O Sequencing_data_lyrata ${accession}
        done < inputList.txt

`inputList.txt` would look like this:

   SRR2040789    
   SRR2040805    
   SRR2040804    


