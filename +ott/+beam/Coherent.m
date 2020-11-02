classdef Coherent < ott.beam.Array
% Coherent array of beams.
% Inherits from :class:`Array`.
%
% Coherent arrays can only contain beams and coherent arrays of beams.
% Adding a incoherent beam to the array raises an error.
%
% Properties
%   - data        -- Array of coherent beams

  methods
    function beam = Coherent(data)
      % Construct a coherent array of beams
      %
      % Usage
      %   beam = Coherent([beam1, beam2, beam3, ...])
      %
      %   beam = Coherent(incoherent_array)
      %   Converts a incoherent array to a coherent array.

      if numel(data) == 1 && isa(data, 'ott.beam.Incoherent')
        data = data.data;
      end

      beam.data = data;
    end
  end

  methods (Hidden)
    function validateArray(beam, data)
      % Check for incoherent arrays
      assert(~any(reshape([data.containsIncoherent], [])),
          'ott:beam:Coherent:incoherent_in_coherent', ...
          'Coherent beam arrays must not contain incoherent arrays');
    end
  end
end
