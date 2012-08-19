% homework.sty 
% For formatting homework  
% --- Class structure: identification part 
% --- 
\ProvidesClass{homework}[9/17/2006 version 1.00 Homework] 
\NeedsTeXFormat{LaTeX2e} 

% --- Class structure: declaration of options part 
% --- 
% This class extends the article class 
% Read all the documentclass options; pass them to article, 
% unless the file "<currentoption>.min" exists, then it is loaded 
\DeclareOption*{\InputIfFileExists{\CurrentOption.hwk}{}{% 
\PassOptionsToClass{\CurrentOption}{article}}} 


% --- Class structure: execution of options part 
% --- 
\ProcessOptions \relax 
% --- Class structure: declaration of options part 
% --- 
\LoadClass{article} 

% Page layout 
\RequirePackage[left=1in,right=1in,top=1in,bottom=1in]{geometry} 


\usepackage{fancyhdr}           % For LaTeX2.09 you should specify {fancyhdr} 
                                % in the \documentstyle instead of the 
                                % \usepackage command
\usepackage{ifthen}

% headers for all pages
\pagestyle{fancy}
\fancyhf{}                      % clear all header and footer fields
\renewcommand{\headrulewidth}{1pt}   
\renewcommand{\footrulewidth}{1pt}   

% if short title (#1) does not exits, use long title (#2) for header, else use #1
\newcommand\settitle[2][]{%
 \title{#2}
 \ifthenelse{\equal{#1}{}}%
  {
  	\fancyhead[C]{#2}
  	\fancyhead[R]{\thepage}
  }%
  {
  	\fancyhead[C]{#1}
  	\fancyhead[R]{\thepage}
  }%
}

% use short author list (#1) for header
\newcommand\setauthors[2]{%
 \author{#2}
  {\fancyhead[L]{#1}}%
}

% header for first page
\fancypagestyle{plain}{%
\fancyhf{}  
% URL suffix will be filled in by Dagstuhl             
%\fancyfoot[L]{\scriptsize Dagstuhl Seminar Proceedings <seminar-number>\\ http://drops.dagstuhl.de/opus/volltexte/2005/}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}}