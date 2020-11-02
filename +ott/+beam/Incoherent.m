class Incoherent < ott.beam.Array
% Incoherent array of beams.
% Inherits from :class:`Array`.
%
% Properties
%   - data        -- Array of incoherent beams

  methods
    function beam = Incoherent(data)
      % Construct a incoherent array of beams
      %
      % Usage
      %   beam = Incoherent([beam1, beam2, ...])
      %
      %   beam = Incoherent(coherent_array)
      %   Converts a coherent array to an incoherent array.

      if numel(data) == 1 && isa(data, 'ott.beam.Coherent')
        data = data.data;
      end

      beam.data = data;
    end
  end

  methods (Hidden)
    function validateArray(~, ~)
      % Nothing to do
    end
  end
end

