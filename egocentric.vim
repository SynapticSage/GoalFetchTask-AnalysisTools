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
let Lf_PopupColorscheme = "gruvbox_material"
let UltiSnipsJumpForwardTrigger = "<C-j>"
let VM_Extend_hl = "Visual"
let Tlist_Ctags_Cmd = "ctags"
let DevIconsEnableDistro =  1 
let Taboo_tabs = "1	Notes\n2	Main\n"
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
let Bold = ""
let Tlist_Exit_OnlyWindow =  0 
let Tlist_Max_Tag_Length =  10 
let UltiSnipsListSnippets = "<c-tab>"
let Tlist_Auto_Open =  0 
let DevIconsEnableFolderExtensionPatternMatching =  0 
let Italic = ""
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
$argadd analyses.m
$argadd \+coding/poissonModel_given_tuningCurveModel.m
$argadd \+coding/+field/combine.m
$argadd \+coding/+field/place.m
$argadd \+coding/+field/calc.m
$argadd \+coding/+test/vectorCoding.m
$argadd \+coding/tuningModel.m
$argadd \+coding/+control/nullFiring.m
$argadd \+coding/+vectorCells/explainHeadDir.m
$argadd \+coding/+vectorCells/+jercog/gpumodel.m
$argadd \+coding/+vectorCells/+jercog/optmodel.m
$argadd \+coding/+vectorCells/+jercog/optgpumodel.m
$argadd \+coding/+vectorCells/+jercog/+table/modelreal.m
$argadd \+coding/+vectorCells/+jercog/+table/+field/nonmarginal.m
$argadd \+coding/+vectorCells/+jercog/+table/+field/standard.m
$argadd \+coding/+vectorCells/+jercog/+table/findmarginals.m
$argadd \+coding/+vectorCells/+jercog/+table/params.m
$argadd \+coding/+vectorCells/+jercog/model.m
$argadd \+coding/+vectorCells/+jercog/allNeurons_sparse.m
$argadd \+coding/+vectorCells/+jercog/dense.m
$argadd \+coding/+vectorCells/+jercog/allNeurons.m
$argadd \+coding/+vectorCells/+jercog/allNeurons_sparse_lean.m
$argadd \+coding/+vectorCells/+jercog/+plot/referenceDistribution.m
$argadd \+coding/+vectorCells/explainPlace.m
$argadd \+coding/+vectorCells/+sarel/rayleigh.m
$argadd \+coding/+vectorCells/+sarel/occupancyNormalize.m
$argadd \+coding/+vectorCells/+sarel/sarel.m
$argadd \+coding/+vectorCells/+sarel/table.m
$argadd \+coding/+vectorCells/+sarel/binning.m
$argadd \+coding/+vectorCells/+sarel/+table/modelreal.m
$argadd \+coding/+vectorCells/+sarel/+table/findmarginals.m
$argadd \+coding/+vectorCells/+sarel/+table/params.m
$argadd \+coding/+vectorCells/+sarel/+table/standardfields.m
$argadd \+coding/+vectorCells/+sarel/tuning.m
$argadd \+coding/+vectorCells/+sarel/rayleighCompute.m
$argadd \+coding/+table/info.m
$argadd \+coding/+table/summarize.m
$argadd \+coding/+plots/jercog.m
$argadd \+coding/+futurepast/+dotson/fieldMI.m
$argadd \+dotson/yotsen.m
$argadd \+coding/+futurepast/+dotson/planofattack.txt
$argadd \+dotson/combinefields.m
$argadd \+coding/+futurepast/+dotson/fieldentropy.m
$argadd \+dotson/calcfield.m
$argadd \+coding/+futurepast/+dotson/+plot/overallDistributions.m
$argadd \+coding/+futurepast/+dotson/+plot/cellDistributions.m
$argadd \+coding/+file/load.m
$argadd \+coding/+file/exist.m
$argadd \+coding/+file/+combineMethods/jercog.m
$argadd \+coding/+file/+combineMethods/sarel.m
$argadd \+coding/+file/+combineMethods/filename.m
$argadd \+coding/+file/filename.m
$argadd \+coding/+file/save.m
$argadd \+units/+clean/+sparse/addMissingNeurons.m
$argadd \+units/getRateMatrix.m
$argadd \+units/getMarkMatrix.m
$argadd \+units/atBehavior_singleCell.m
$argadd \+units/atBehavior.m
$argadd \+gbehavior/+clean/removeduplicatetimes.m
$argadd \+gbehavior/lookup.m
$argadd \+gbehavior/mazebounds.m
$argadd \+tabfilt/behavior.m
$argadd \+tabfilt/spikesbeh.m
$argadd \+tabfilt/spikes.m
$argadd egocentric.vim
$argadd \+maze/imagesc.m
tabnew
tabrewind
edit \+units/sparseToDense.m
argglobal
if bufexists("\+units/sparseToDense.m") | buffer \+units/sparseToDense.m | else | edit \+units/sparseToDense.m | endif
if &buftype ==# 'terminal'
  silent file \+units/sparseToDense.m
endif
balt ~/Documents/Notes/Overarch
let s:l = 211 - ((38 * winheight(0) + 20) / 41)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 211
normal! 0
lcd ~/Code/analysis
tabnext
edit ~/Code/analysis/analyses.m
let s:save_splitbelow = &splitbelow
let s:save_splitright = &splitright
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
let &splitbelow = s:save_splitbelow
let &splitright = s:save_splitright
wincmd t
let s:save_winminheight = &winminheight
let s:save_winminwidth = &winminwidth
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 61 + 61) / 122)
exe 'vert 2resize ' . ((&columns * 60 + 61) / 122)
argglobal
balt ~/Code/analysis/+units/atBehavior.m
let s:l = 1 - ((0 * winheight(0) + 18) / 37)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
lcd ~/Code/analysis
wincmd w
argglobal
if bufexists("~/Code/analysis/+units/atBehavior_singleCell.m") | buffer ~/Code/analysis/+units/atBehavior_singleCell.m | else | edit ~/Code/analysis/+units/atBehavior_singleCell.m | endif
if &buftype ==# 'terminal'
  silent file ~/Code/analysis/+units/atBehavior_singleCell.m
endif
balt ~/Code/analysis/+units/atBehavior.m
let s:l = 34 - ((11 * winheight(0) + 18) / 37)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 34
normal! 016|
lcd ~/Code/analysis
wincmd w
exe 'vert 1resize ' . ((&columns * 61 + 61) / 122)
exe 'vert 2resize ' . ((&columns * 60 + 61) / 122)
tabnext 2
badd +1 ~/Documents/Notes/Overarch
badd +1 ~/Code/analysis/analyses.m
badd +1 ~/Code/analysis/+coding/poissonModel_given_tuningCurveModel.m
badd +1 ~/Code/analysis/+coding/+field/combine.m
badd +1 ~/Code/analysis/+coding/+field/place.m
badd +132 ~/Code/analysis/+coding/+field/calc.m
badd +1 ~/Code/analysis/+coding/+test/vectorCoding.m
badd +1 ~/Code/analysis/+coding/+control/nullFiring.m
badd +1 \+coding/+vectorCells/explainHeadDir.m
badd +1 \+coding/+vectorCells/+jercog/gpumodel.m
badd +1 \+coding/+vectorCells/+jercog/optmodel.m
badd +1 \+coding/+vectorCells/+jercog/optgpumodel.m
badd +1 \+coding/+vectorCells/+jercog/+table/modelreal.m
badd +1 \+coding/+vectorCells/+jercog/+table/+field/nonmarginal.m
badd +1 \+coding/+vectorCells/+jercog/+table/+field/standard.m
badd +1 \+coding/+vectorCells/+jercog/+table/findmarginals.m
badd +1 \+coding/+vectorCells/+jercog/+table/params.m
badd +1 \+coding/+vectorCells/+jercog/model.m
badd +1 \+coding/+vectorCells/+jercog/allNeurons_sparse.m
badd +1 \+coding/+vectorCells/+jercog/dense.m
badd +1 \+coding/+vectorCells/+jercog/allNeurons.m
badd +1 \+coding/+vectorCells/+jercog/allNeurons_sparse_lean.m
badd +1 \+coding/+vectorCells/+jercog/+plot/referenceDistribution.m
badd +1 \+coding/+vectorCells/explainPlace.m
badd +1 \+coding/+vectorCells/+sarel/rayleigh.m
badd +1 \+coding/+vectorCells/+sarel/occupancyNormalize.m
badd +1 \+coding/+vectorCells/+sarel/sarel.m
badd +1 \+coding/+vectorCells/+sarel/table.m
badd +1 \+coding/+vectorCells/+sarel/binning.m
badd +1 \+coding/+vectorCells/+sarel/+table/modelreal.m
badd +1 \+coding/+vectorCells/+sarel/+table/findmarginals.m
badd +1 \+coding/+vectorCells/+sarel/+table/params.m
badd +1 \+coding/+vectorCells/+sarel/+table/standardfields.m
badd +1 \+coding/+vectorCells/+sarel/tuning.m
badd +1 \+coding/+vectorCells/+sarel/rayleighCompute.m
badd +1 ~/Code/analysis/+coding/+table/info.m
badd +1 ~/Code/analysis/+coding/+table/summarize.m
badd +1 ~/Code/analysis/+coding/+plots/jercog.m
badd +1 \+coding/+futurepast/+dotson/fieldMI.m
badd +1 \+dotson/yotsen.m
badd +1 \+coding/+futurepast/+dotson/planofattack.txt
badd +1 \+coding/+futurepast/+dotson/fieldentropy.m
badd +1 \+coding/+futurepast/+dotson/+plot/overallDistributions.m
badd +1 \+coding/+futurepast/+dotson/+plot/cellDistributions.m
badd +6 ~/Code/analysis/+coding/+file/load.m
badd +1 ~/Code/analysis/+coding/+file/exist.m
badd +1 ~/Code/analysis/+coding/+file/+combineMethods/jercog.m
badd +1 ~/Code/analysis/+coding/+file/+combineMethods/sarel.m
badd +1 ~/Code/analysis/+coding/+file/+combineMethods/filename.m
badd +20 ~/Code/analysis/+coding/+file/filename.m
badd +43 ~/Code/analysis/+coding/+file/save.m
badd +32 ~/Code/analysis/+units/+clean/+sparse/addMissingNeurons.m
badd +10 ~/Code/analysis/+units/getRateMatrix.m
badd +1 ~/Code/analysis/+units/getMarkMatrix.m
badd +4 ~/Code/analysis/+units/atBehavior_singleCell.m
badd +30 ~/Code/analysis/+units/atBehavior.m
badd +1 ~/Code/analysis/+gbehavior/+clean/removeduplicatetimes.m
badd +124 ~/Code/analysis/+gbehavior/lookup.m
badd +1 ~/Code/analysis/+gbehavior/mazebounds.m
badd +1 ~/Code/analysis/+tabfilt/behavior.m
badd +1 ~/Code/analysis/+tabfilt/spikesbeh.m
badd +1 ~/Code/analysis/+tabfilt/spikes.m
badd +1 ~/Code/analysis/egocentric.vim
badd +1 ~/Code/analysis/+maze/imagesc.m
badd +51 ~/Code/analysis/+units/+shuffle/conditional_time.m
badd +5 ~/Code/analysis/tmp.m
badd +1 ~/Code/analysis/+connectivity/+CoNNECT/modules/makeCC.py
badd +1 ~/Code/analysis/+coding/+sarel/+table/tuning.m
badd +1 ~/Code/analysis/+coding/+sarel/computeGoalPlaceIndex.m
badd +31 ~/Code/analysis/+coding/+sarel/binning.m
badd +1 ~/Code/analysis/+connectivity/+CoNNECT/modules/sources/crosscorrelogram.pyx
badd +38 ~/Code/analysis/+connectivity/+CoNNECT/modules/crosscorrelogram_gpu.py
badd +1 ~/anaconda3/lib/python3.9/site-packages/numba/core/decorators.py
badd +1 ~/Code/analysis/+units/+shuffle/yield.m
badd +1 ~/Code/analysis/+coding/+jercog/Plot.m
badd +1 ~/Code/analysis/+coding/+sarel/Plot.m
badd +1 ~/Code/analysis/plotanalyses.m
badd +1 ~/Code/analysis/+coding/+sarel/+helper/workflow.m
badd +63 ~/Code/analysis/+connectivity/+CoNNECT/estimate.py
badd +1 ~/Code/analysis/+gbehavior/+phy/makeNPYfiles.m
badd +5 ~/Code/analysis/+gbehavior/+phy/vars2d.m
badd +1 ~/Code/analysis/+gbehavior/+phy/atomicNPYfile.m
badd +1 ~/translate.txt
badd +1 ~/Code/projects/pointprocess/replay_trajectory_classification/replay_trajectory_classification/multiunit_likelihood_gpu.py
badd +1 ~/Code/analysis/+connectivity/+CoNNECT/modules/utils.py
badd +1 ~/Code/analysis/+connectivity/+glmcc/GLMCC/glmcc_fitting.py
badd +2 ~/.vim/plugged/matlab.vim/doc/matlab.txt
badd +1 ~/Code/analysis/+coding/+sarel/computeMaxMeanIndices.m
badd +1 ~/Code/analysis/+coding/+sarel/sarel.m
badd +1 ~/.vim/plugged/vim-notes/misc/notes/user/jercog
badd +1 ~/Code/analysis
badd +57 ~/Code/analysis/+coding/+jercog/+table/params.m
badd +1 ~/Code/Src_Matlab/CircStat2012a/circ_mean.m
badd +1 ~/Code/analysis/+coding/+sarel/+helper/GPUworkflow.m
badd +1 ~/Code/Src_Matlab/CircStat2012a/circ_r.m
badd +1 ~/Code/analysis/jercog.m
badd +1 ~/Code/analysis/+coding/+jercog/+plot/referenceDistribution.m
badd +1 ~/Code/analysis/+coding/+vectorCells/explainHeadDir.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/gpumodel.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/optmodel.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/optgpumodel.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/+table/modelreal.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/+table/+field/nonmarginal.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/+table/+field/standard.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/+table/findmarginals.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/+table/params.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/model.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/allNeurons_sparse.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/dense.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/allNeurons.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/allNeurons_sparse_lean.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+jercog/+plot/referenceDistribution.m
badd +1 ~/Code/analysis/+coding/+vectorCells/explainPlace.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/rayleigh.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/occupancyNormalize.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/sarel.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/table.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/binning.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/+table/modelreal.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/+table/findmarginals.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/+table/params.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/+table/standardfields.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/tuning.m
badd +1 ~/Code/analysis/+coding/+vectorCells/+sarel/rayleighCompute.m
badd +1 ~/Documents/Notes/goalmaze.sarel
badd +1 ~/Code/analysis/+coding/+jercog/allNeurons_sparse_lean.m
badd +1 ~/Code/pipeline/ry_pipeline/+lfpLib/+noise/coherence.m
badd +1 man://ip(8)
badd +1 ~/.vim/plugged/DrawIt/doc/DrawIt.txt
badd +1 ~/Code/analysis/+coding/+vectorCells/sarel.m
badd +1 /Volumes/MATLAB-Drive/Shared/+table/+behavior/lookup.m
badd +1 ~/Code/pipeline/ry_pipeline/+rawLib/+mod/addmarks.m
badd +1 \+file/+combineMethods/sarel.m
badd +1 /Volumes/MATLAB-Drive/linearized/+tidyData/fromNdb.m
badd +1 man://end(3)
badd +18 /Volumes/MATLAB-Drive/linearized/+util/+job/waitQueueSize.m
badd +4 /Volumes/MATLAB-Drive/linearized/+util/+job/continueIfExist.m
badd +1 man://utime(2)
badd +34 ~/Code/Src_Matlab/CircStat2012a/circ_rtest.m
badd +1 man://cluster(1)
badd +1 ~/.vim/plugged/vim-notes/misc/notes/user/LEAN
badd +1 ~/.vim/plugged/vim-notes/misc/notes/user/NOTE
badd +1 ~/.vim/plugged/vim-notes/misc/notes/user/New\ note
badd +1 ~/.vim/plugged/vim-notes/misc/notes/user/allNeurons_sparse_LEAN
badd +1 ~/.vim/plugged/vim-notes/misc/notes/user/allNeurons_sparse_LEA
badd +1 ~/.vim/plugged/vim-notes/misc/notes/user/jercog_replication
badd +1 ~/Code/analysis/Tagged\ Notes
badd +1 ~/Code/analysis/-type
badd +1 ~/Code/analysis/f\'
badd +1 ~/Code/analysis/+coding
badd +1 ~/Code/analysis/f
badd +1 ~/Code/analysis/-name
badd +1 ~/Code/analysis/\`find\ .\ -type\ f\ -name
badd +3 ~/Code/analysis/+coding/+jercog/+table/findmarginals.m
badd +6 ~/Code/analysis/+coding/+jercog/+table/+field/standard.m
badd +1 ~/Code/analysis/+coding/+jercog/+table/+field/nonmarginal.m
badd +3 ~/Code/analysis/+coding/+sarel/+table/+field/standard.m
badd +55 ~/Code/analysis/+coding/+sarel/tuning.m
badd +1 ~/Code/analysis/+coding/+sarel/table.m
badd +1 ~/Code/analysis/+coding/+file/matfile.m
badd +15 ~/Code/analysis/+coding/+sarel/+table/findmarginals.m
badd +4 ~/Code/analysis/+coding/+sarel/+table/+field/standardBin.m
badd +57 ~/Code/analysis/+coding/+field/scaffold.m
badd +18 ~/Code/analysis/+coding/+field/scaffold2binning.m
badd +26 ~/Code/analysis/+coding/+sarel/+get/goal.m
badd +4 ~/Code/analysis/+coding/+sarel/+get/place.m
badd +1 ~/Code/analysis/+coding/+sarel/computeDirectionalityIndex.m
badd +1 ~/Code/analysis/+units/+shuffle/time.m
badd +2 ~/Code/analysis/+units/+shuffle/unit.m
badd +1 ~/Code/analysis/+units/+shuffle/condition.m
badd +1 ~/Code/analysis/+units/+shuffle/+helper/preallocateBehaviorShifts.m
badd +3 ~/Code/analysis/+coding/+file/shufflematfilename.m
badd +3 ~/Code/analysis/+coding/+file/shufflematfile.m
badd +147 ~/Code/analysis/+connectivity/+glmcc/GLMCC/Est_Data.py
badd +22 ~/Code/analysis/+connectivity/+CoNNECT/fitting.m
badd +43 ~/Code/analysis/+connectivity/makecelltmpfiles.m
badd +4963 ~/Code/analysis/+connectivity/+CoNNECT/modules/sources/crosscorrelogram.c
badd +5 ~/.vimpyter_views/Numba_GPU_KDE.ipynb
badd +1 ~/Code/projects/pointprocess/replay_trajectory_classification/replay_trajectory_classification/misc_cuda2.py
badd +1 ~/Code/analysis/__doc__
badd +1 ~/Code/Src_Matlab/ry_HomeBrew_Utility/datadef.m
badd +50 ~/Code/pipeline/preprocess/preprocess_RY16.m
badd +1 ~/Code/analysis/+gbehavior/+phy/vars1d.m
badd +1 ~/Code/pipeline/ry_ms4_pipeline/convert_ml_to_FF_withMultiunitAndMarks.m
badd +197 ~/Code/pipeline/ry_ms4_pipeline/ml_process_animal.m
badd +90 ~/Code/pipeline/ry_ms4_pipeline/ml_sort_on_segs.m
badd +27 ~/Code/Src_Matlab/ry_Utility/npy-matlab/npy-matlab/constructNPYheader.m
badd +1 ~/Code/analysis/+units/+shuffle/conditional_time_yield.m
badd +20 ~/Code/analysis/+coding/+futurepast/main.m
badd +77 ~/Code/analysis/+coding/+futurepast/fieldMI.m
badd +1 ~/Documents/Notes/Shuffle\ Change
badd +1 ~/Documents/Notes/New\ note
badd +10 ~/Code/analysis/+coding/+file/shuffleparquetfolder.m
badd +18 ~/Code/analysis/+units/+shuffle/get.m
badd +10 ~/Code/analysis/+coding/+file/+shuffle/ParquetToMatfile.m
badd +1 ~/Code/analysis/+units/+shuffle/+helper/measuretimeperiods.m
badd +1 ~/Code/analysis/+units/+shuffle/+helper/calculateshifts.m
badd +1 ~/Code/analysis/+units/+shuffle/+helper/prepareOuts.m
badd +67 ~/Code/analysis/+units/+shuffle/+helper/behaviorbasedshuffle.m
badd +1 ~/Code/analysis/+units/+shuffle/+helper/cache.m
badd +33 /usr/local/MATLAB/R2021b/toolbox/matlab/datatypes/tabular/struct2table.m
badd +2 ~/Code/analysis/tmp_histcount.m
badd +211 ~/Code/analysis/+units/sparseToDense.m
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0 && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToOFIAc
let &winminheight = s:save_winminheight
let &winminwidth = s:save_winminwidth
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
