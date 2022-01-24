function folder = datafolder(varargin)

folder = fullfile(coding.file.projectfolder, 'data', varargin{:});
