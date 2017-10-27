% plot_solution.m, Version 1.0
% (c) Oliver J. Sutton, University of Leicester, 2016.
% Distributed for academic and educational purposes only - all rights of reproduction reserved.
%
% If you find this code helpful, please reference the associated paper:
%			Sutton, O.J., "The virtual element method in 50 lines of Matlab", 2016, Submitted.
% which also contains the documentation of the code.
function plot_solution(mesh, solution)
% Plot the vertex values of the virtual element solution
	title('Approximate Solution');
	maxNumVertices = max(cellfun(@numel, mesh.elements));
	padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
	elements = cellfun(padFunc, mesh.elements, 'UniformOutput', false);
	elements = vertcat(elements{:});
	data = [mesh.vertices, solution];
	patch('Faces', elements, 'Vertices', data,...
							'FaceColor', 'flat', 'CData', solution / max(abs(solution)));
	%axis('square')
	xlim([min(mesh.vertices(:, 1)) - 0.1, max(mesh.vertices(:, 1)) + 0.1])
	ylim([min(mesh.vertices(:, 2)) - 0.1, max(mesh.vertices(:, 2)) + 0.1])
	zlim([min(solution) - 0.1, max(solution) + 0.1])
	xlabel('x'); ylabel('y'); zlabel('u');
end