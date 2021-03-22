#!/bin/bash
usage=$(cat <<-END
    runit [-h] [-d] [-cl] [--dry]

    This utility will run all R scripts in a directory, if numbered
    then in the order they are numbered. Defaults to running all scripts
    in subdirectories of directory R:
        runit -d "R/*/*"
    You can skip a script by including a commented line with the word "skipit" in the script.
    accepted options are:
    -h, --help          show help message and exit
    -d, --dir           the directory path, quoted
    -cl, --cluster      whether to run the cluster scripts, defaults to not running them,
                        if passed, will ask you whether vpn-ed in to della
    --dry               dryrun, just prints the scripts and whether they are run
                        locally or on cluster
    -q, --quiet         whether to show messages, warnings, and errors from R; defaults to showing all
                        pass the arg to suppress these
    --printErrors       whether to print errors regardless of quiet argument
    --force             if you want to force the skipped jobs to run
END
)

chs=$(cat <<-END
How do you want to run the cluster jobs?
    1 | local (in parallel with available cores on your local computer)
    2 | serial (run without a parallel backend locally)
    3 | della (you have access to the della cluster and have set up the sub utility cmd)
    4 | skip  (do not run the cluster jobs.)
Your choice:
END
)

# Set up
quiet=0
dryrun=0
cl=0
printerrors=0
skip=1
while [ "$1" != "" ]; do
    case $1 in
    -d | --dir )           shift
                            FILES=$1
                            ;;
    --dry )                 dryrun=1
                            ;;
    -q | --quiet )          quiet=1
                            ;;
    --printErrors )         printerrors=1
                            ;;
    -cl | --cluster )       cl=1
                            ;;
    --force )               skip=0
                            ;;
    -h | --help )           echo "$usage"
                            exit
                            ;;
    * )                     echo "$usage"
                            exit 1
    esac
    shift
done

if [ "$quiet" = "0" ];
then
    out=/dev/tty
else
    if [ "$printerrors" = "0" ];
    then
        out=/dev/null
    else
        out=$(mktemp /tmp/log.XXXXXX)
        trap "rm -f $out" 0 3 15 6
    fi
fi

# Colors
Red='\033[0;31m'          # Red
Blue='\033[0;34m'         # Blue
NC='\033[0m' # No Color
BCyan='\033[1;36m' # Bold Cyan
BRed='\033[1;31m' # Bold Red

# Processing stdout & stderr

if [ "$cl" = "1" ] && [ "$dryrun" = "0" ];
then
read -p "$chs" choice
    case "$choice" in
        1|local ) arg=local; echo "Running locally and in parallel!";;
        2|serial ) arg=serial;echo "Running serially!";;
        3|della ) arg=della ; echo "Running on della!";;
        4|skip ) cl=0; echo "Skipping cluster jobs!";;
        * ) echo "invalid";exit;;
    esac
fi

for f in $FILES
do
    if grep -q "skipit" "$f" ;
    then
    echo -e "${Red}Skipping job: $f ${NC}"
    continue
    fi

    if grep -q "sub_cmd" "$f";
    then
        if [ "$cl" = "0" ];
        then
            echo "cluster job not run: $f"
        else
            if [ "$arg" = "della"  ];
            then
                cmd=$(grep "sub_cmd"  "$f" | cut -f 2 -d =) # sub args from script
                if [ "$dryrun" = "1" ];
                then
                    echo "cluster cmd: sub $cmd -sp $f"
                else
                    echo -e "${Blue}Cluster job: $f${NC}"
                    sub $cmd -sp $f
                fi
            else
                if [ "$dryrun" = "1" ];
                then
                    echo -e "${Blue}Local job ($arg): $f ${NC}"
                else
                    echo -e "${Blue}Local job ($arg): $f ${NC}"
                    if Rscript --vanilla "$f" "$arg" &> $out;
                    then
                        echo  -e "${BCyan}$f completed.${NC}"
                    else
                        echo -e "${BRed}$f did not complete!${NC}"
                        if [ "$printerrors" = "1" ];
                        then
                            grep "Error" $out --color
                        fi
                    fi
                fi
            fi
        fi
    else
        if [ "$dryrun" = "1" ];
        then
            echo -e "${Blue}Local job: $f ${NC}"
        else
            echo -e "${Blue}Local job: $f ${NC}"
            if Rscript --vanilla "$f" &> $out;
            then
                echo -e "${BCyan}$f completed.${NC}"
            else
                echo -e "${BRed}$f did not complete!${NC}"
                if [ "$printerrors" = "1" ];
                then
                    grep "Error" $out --color
                fi
            fi
        fi
    fi
done

echo "" # for trap even if you CTRL + C your way out of all the jobs
