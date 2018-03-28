fields = fieldnames(flightdata);
T = table
for i = 1:numel(fields)
  a = fields{i};
  a
  b = flightdata.(fields{i}).data;
  if i == 48
      b = b.'
  end
  T.(a) = b;
end% see the new result...
writetable(T, 'flightdata.txt')

