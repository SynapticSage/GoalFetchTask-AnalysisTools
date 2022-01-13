let SessionLoad = 1
let s:so_save = &g:so | let s:siso_save = &g:siso | setg so=0 siso=0 | setl so=-1 siso=-1
let v:this_session=expand("<sfile>:p")
let Tlist_Show_One_File =  0 
let WebDevIconsNerdTreeGitPluginForceVAlign =  1 
let UltiSnipsRemoveSelectModeMappings =  1 
let WebDevIconsUnicodeDecorateFolderNodesDefaultSymbol = ""
let UltiSnipsExpandTrigger = "<C-j>"
let Tlist_GainFocus_On_ToggleOpen =  0 
let Tlist_Auto_Update =  1 
let Lf_StlColorscheme = "one"
let UltiSnipsDebugHost = "localhost"
let Tlist_File_Fold_Auto_Close =  0 
let UltiSnipsEditSplit = "vertical"
let VM_Cursor_hl = "Cursor"
let WebDevIconsTabAirLineBeforeGlyphPadding = " "
let Tlist_Sort_Type = "order"
let Tlist_Enable_Fold_Column =  1 
let Tlist_Compact_Format =  0 
let UltiSnipsJumpForwardTrigger = "<C-j>"
let VM_Extend_hl = "Visual"
let Tlist_Ctags_Cmd = "ctags"
let DevIconsEnableDistro =  1 
let Taboo_tabs = "1	Sarel:Main\n2	Sarel\n"
let WebDevIconsUnicodeDecorateFileNodes =  1 
let UltiSnipsPMDebugBlocking =  0 
let Tlist_Display_Tag_Scope =  1 
let Tlist_Highlight_Tag_On_BufEnter =  1 
let Tlist_Max_Submenu_Items =  20 
let Tlist_WinWidth =  30 
let Tlist_Close_On_Select =  0 
let NERDTreeUpdateOnCursorHold =  1 
let DevIconsDefaultFolderOpenSymbol = ""
let UltiSnipsJumpBackwardTrigger = "<C-k>"
let UltiSnipsDebugPort =  8080 
let Tlist_Show_Menu =  0 
let Tlist_Use_SingleClick =  0 
let UltiSnipsUsePythonVersion =  3 
let UltiSnipsEnableSnipMate =  1 
let DevIconsEnableFoldersOpenClose =  0 
let DevIconsEnableFolderPatternMatching =  1 
let VM_Insert_hl = "Cursor"
let TagList_title = "__Tag_List__"
let WebDevIconsTabAirLineAfterGlyphPadding = ""
let WebDevIconsUnicodeDecorateFileNodesDefaultSymbol = ""
let NERDTreeGitStatusUpdateOnCursorHold =  1 
let WebDevIconsUnicodeDecorateFolderNodesExactMatches =  1 
let WebDevIconsUnicodeByteOrderMarkerDefaultSymbol = ""
let DevIconsAppendArtifactFix =  0 
let WebDevIconsUnicodeDecorateFolderNodes =  1 
let UltiSnipsDebugServerEnable =  0 
let WebDevIconsUnicodeGlyphDoubleWidth =  1 
let Tlist_Exit_OnlyWindow =  0 
let Tlist_Max_Tag_Length =  10 
let UltiSnipsListSnippets = "<c-tab>"
let Tlist_Auto_Open =  0 
let DevIconsEnableFolderExtensionPatternMatching =  0 
let Tlist_Auto_Highlight_Tag =  1 
let Tlist_Use_Horiz_Window =  0 
let WebDevIconsNerdTreeAfterGlyphPadding = " "
let Tlist_Process_File_Always =  0 
let ScreenImpl = "Tmux"
let VM_Mono_hl = "Cursor"
let DevIconsArtifactFixChar = " "
let Tlist_WinHeight =  10 
let Tlist_Inc_Winwidth =  1 
let WebDevIconsNerdTreeBeforeGlyphPadding = " "
let WebDevIconsUnicodeDecorateFolderNodesSymlinkSymbol = ""
let Tlist_Display_Prototype =  0 
let Tlist_Use_Right_Window =  0 
silent only
silent tabonly
cd ~/Code/analysis
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
argglobal
%argdel
$argadd oh
tabnew
tabnew
tabrewind
edit analyses.m
argglobal
balt \+units/atBehavior_singleCell.m
let s:l = 91 - ((0 * winheight(0) + 20) / 40)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 91
normal! 0
lcd ~/Code/analysis
tabnext
edit ~/Code/analysis/+coding/+sarel/computeContrastIndex.m
argglobal
balt ~/Code/analysis/+coding/+sarel/computeMaxMeanIndices.m
let s:l = 19 - ((18 * winheight(0) + 21) / 43)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 19
normal! 026|
lcd ~/Code/analysis
tabnext
edit ~/Code/analysis/+coding/+field/calc.m
argglobal
balt ~/Code/analysis/+coding/+futurepast/main.m
let s:l = 1 - ((0 * winheight(0) + 17) / 35)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
lcd ~/Code/analysis
tabnext 3
badd +8 ~/Code/analysis/+coding/+sarel/computeDirectionalityIndex.m
badd +0 ~/Code/analysis/oh
badd +128 ~/Code/analysis/analyses.m
badd +0 ~/Code/analysis/+coding/+futurepast/shuffcorrect.asv
badd +1 ~/Code/analysis/+units/+shuffle/get.m
badd +44 ~/Code/analysis/+units/+shuffle/+helper/behaviorbasedshuffle.m
badd +1 ~/Code/analysis/+coding/+file/+combineMethods/sarel.m
badd +1 ~/Code/analysis/+coding/+sarel/main.m
badd +1 ~/Code/analysis/+coding/+sarel/computeGoalPlaceIndex.m
badd +40 ~/Code/analysis/+coding/+futurepast/main.m
badd +97 ~/Code/analysis/+coding/+field/calc.m
badd +1 ~/Code/analysis/+coding/+field/scaffold.m
badd +1 ~/Code/analysis/+coding/+sarel/binning.m
badd +1 ~/Code/analysis/+coding/+field/scaffold2binning.m
badd +8 ~/Code/analysis/+coding/+sarel/+helper/workflow.m
badd +25 ~/Code/analysis/+coding/+sarel/tuning.m
badd +0 ~/Code/analysis/+coding/+sarel/Plot.m
badd +7 ~/Code/analysis/+units/indexToBehavior.m
badd +2 ~/Code/analysis/+coding/+sarel/+get/goal.m
badd +44 ~/Code/analysis/+coding/+sarel/+get/vonmises.m
badd +1 ~/Code/analysis/+coding/+sarel/computeFitVonMises.m
badd +1 ~/Code/analysis/+coding/+sarel/computeMaxMeanIndices.m
badd +29 ~/Code/Src_Matlab/ry_Utility/jpca-extension/bin/circ_mean.m
badd +24 ~/Code/analysis/+units/atBehavior_singleCell.m
badd +1 ~/Code/analysis/+coding/+sarel/+get/place.m
badd +22 ~/Code/analysis/+coding/+sarel/computeContrastIndex.m
badd +60 ~/Code/analysis/+coding/+sarel/analysis.m
badd +72 ~/Code/analysis/+coding/+futurepast/analysis.m
badd +2 ~/Code/analysis/+coding/+jercog/analysis.m
badd +21 ~/Code/analysis/+coding/+sarel/+get/contrast.m
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0 && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToOFAIc
let s:sx = expand("<sfile>:p:r")."x.vim"
if filereadable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &g:so = s:so_save | let &g:siso = s:siso_save
set hlsearch
let g:this_session = v:this_session
let g:this_obsession = v:this_session
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
