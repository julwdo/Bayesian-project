\NeedsTeXFormat{LaTeX2e}
\IfFileExists{Class/julia}{\ProvidesClass{Class/julia}}
    \LoadClass{article}
    \typeout{*************}
    \typeout{* julia *}
    \typeout{*************}

% PATHCLASS definition

%% Useful packages
%    \RequirePackage{verbatim}
%    \RequirePackage{etoolbox}
%    \RequirePackage{ifthen}
%    \RequirePackage{enumitem}
%    \RequirePackage{calc}
%
%% Hyperref and url
%    \RequirePackage{url}
%    \RequirePackage{hyperref}
%
%% Font
%    \RequirePackage[T1]{fontenc}
%    \RequirePackage[latin2,latin9,utf8]{inputenc}
%    \RequirePackage{ae}
%
%% Language
%    %
%    \RequirePackage[english]{babel}
%
%% Layout
%    \RequirePackage{geometry}
%    \geometry{paper=a4paper,top=1.5cm,bottom=1.5cm,left=1.5cm,right=1.5cm}
%    \setlength{\headheight}{2pt}
%    \RequirePackage{setspace}
%    \renewcommand{\baselinestretch}{1.15}

% Bibliography
     \RequirePackage[authoryear]{natbib}
    \newcommand{\printbibliography}[1]{
        \nocite{*}
        \bibliographystyle{Class/apa-ds}
        \bibliography{#1}
    }

%% Math
%    \RequirePackage{amsmath}
%    \RequirePackage{amssymb}
%    \RequirePackage{siunitx}
%
%% Tables
%    \RequirePackage{tabularx}
%    \RequirePackage{tabulary}
%    \RequirePackage{multirow}
%    \RequirePackage{rotating}
%    \RequirePackage{float}
%    \RequirePackage{booktabs}
%
%% Graphics
%    \RequirePackage{graphicx}
%    \RequirePackage{color}
%
%% Backup packages
%    %\RequirePackage{cases}
%    %\RequirePackage{pdflscape}
%    %\RequirePackage{afterpage}
%    %\RequirePackage{caption}
%    %\RequirePackage{subcaption}
%    %\RequirePackage{txfonts}
%    %\RequirePackage{mathdots}
%    %\RequirePackage[classicReIm]{kpfonts}
%
%% Tokens
%
%% Running footer and headers
%
%% Formatting abstract
%
%% Formatting authors
%
%% TikZ
%    \RequirePackage{tikz}
%    \usetikzlibrary{calc}
%% Logos and copyright footnotes

\endinput