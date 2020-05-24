function verifyProperties(testCase, meta, newbeam, oldbeam)
% Verifies all properties of new-beam match old-beam
%
% Usage
%   verifyProperties(testCase, meta, newbeam, oldbeam)

  props = {meta.PropertyList(~[meta.PropertyList.Dependent]).Name};
  for p = props
    testCase.verifyEqual(newbeam.(p{1}), oldbeam.(p{1}), p{1});
  end
end
