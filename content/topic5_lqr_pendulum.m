clear; close all; clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
figDir = fullfile(projectRoot, 'figures', 'topic5');
buildDir = fullfile(projectRoot, 'build');
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(buildDir, 'dir'), mkdir(buildDir); end

M = 0.5;
m = 0.2;
b = 0.1;
l = 0.3;
I = 0.006;
g = 9.8;
p = I * (M + m) + M * m * l^2;

A = [0, 1, 0, 0;
     0, -((I + m * l^2) * b) / p, (m^2 * g * l^2) / p, 0;
     0, 0, 0, 1;
     0, -(m * l * b) / p, (m * g * l * (M + m)) / p, 0];
B = [0; (I + m * l^2) / p; 0; m * l / p];
C = [1, 0, 0, 0;
     0, 0, 1, 0];
D = [0; 0];

Q = diag([1000, 1, 100, 1]);
R = 1;
K_lqr = lqr(A, B, Q, R);
Acl = A - B * K_lqr;
Cref = [1, 0, 0, 0];
Nbar = -1 / (Cref * (Acl \ B));

reference = 0.2;
t = linspace(0, 5, 2501)';
sysClosed = ss(Acl, B * Nbar, C, D);
[y, t, x] = lsim(sysClosed, reference * ones(size(t)), t);
cartPosition = y(:, 1);
pendulumAngle = y(:, 2);
controlInput = -x * K_lqr' + Nbar * reference;

riseTimeX = rise_time(t, cartPosition, 0.1 * reference, 0.9 * reference);
settlingTimeX = settling_time(t, cartPosition, reference, 0.02 * abs(reference));
settlingTimePhi = settling_time(t, pendulumAngle, 0, 0.02);
overshootPhiRad = max(abs(pendulumAngle));
closedLoopEigenvalues = eig(Acl);

assert(settlingTimeX < 2, 'Topic 5 RED/GREEN failed: x settling time is not below 2 s.');
assert(settlingTimePhi < 2, 'Topic 5 RED/GREEN failed: phi settling time is not below 2 s.');
assert(riseTimeX < 0.5, 'Topic 5 RED/GREEN failed: x rise time is not below 0.5 s.');
assert(overshootPhiRad < 0.35, 'Topic 5 RED/GREEN failed: phi overshoot is not below 0.35 rad.');

configure_report_graphics();
exportDpi = 300;
plotSize = [100, 100, 1300, 780];
tallPlotSize = [100, 100, 1300, 860];
polePlotSize = [100, 100, 1100, 820];

fig = figure('Color', 'w', 'Units', 'pixels', 'Position', tallPlotSize);
subplot(2, 1, 1);
plot(t, cartPosition, 'Color', [0.00, 0.22, 0.80], 'LineWidth', 2.4); hold on;
yline(reference, 'k--', 'LineWidth', 1.9);
grid on;
xlabel('Time (s)'); ylabel('x (m)');
title('Cart position response');
subplot(2, 1, 2);
plot(t, pendulumAngle, 'Color', [0.85, 0.12, 0.00], 'LineWidth', 2.4);
grid on;
xlabel('Time (s)'); ylabel('Phi (rad)');
title('Pendulum angle response');
export_report_figure(fig, fullfile(figDir, 'state_responses.png'), exportDpi);
close(fig);

fig = figure('Color', 'w', 'Units', 'pixels', 'Position', plotSize);
plot(t, controlInput, 'Color', [0.55, 0.00, 0.55], 'LineWidth', 2.4);
grid on;
xlabel('Time (s)'); ylabel('F (N)');
title('LQR control input');
export_report_figure(fig, fullfile(figDir, 'control_input.png'), exportDpi);
close(fig);

fig = figure('Color', 'w', 'Units', 'pixels', 'Position', polePlotSize);
plot(real(closedLoopEigenvalues), imag(closedLoopEigenvalues), 'x', 'Color', [0.00, 0.22, 0.80], 'LineWidth', 2.8, 'MarkerSize', 12); hold on;
xline(0, 'k-', 'LineWidth', 1.2);
yline(0, 'k-', 'LineWidth', 1.2);
xlim([-9, 1]);
ylim([-6, 6]);
grid on;
sgrid;
lighten_sgrid(gca);
plot(real(closedLoopEigenvalues), imag(closedLoopEigenvalues), 'x', 'Color', [0.00, 0.22, 0.80], 'LineWidth', 2.8, 'MarkerSize', 12);
xlabel('Real axis'); ylabel('Imaginary axis');
title('Closed-loop eigenvalues with damping grid');
export_pole_figure(fig, fullfile(figDir, 'closed_loop_eigenvalues.png'), exportDpi);
close(fig);

create_topic5_simulink_model(projectRoot, A, B, K_lqr, Nbar);

metricsPath = fullfile(buildDir, 'topic5_metrics.txt');
fid = fopen(metricsPath, 'w', 'n', 'UTF-8');
fprintf(fid, 'A =\n'); fprintf(fid, '  %.8f %.8f %.8f %.8f\n', A');
fprintf(fid, 'B = [%.8f %.8f %.8f %.8f]^T\n', B);
fprintf(fid, 'Q = diag([1000, 1, 100, 1])\n');
fprintf(fid, 'R = %.8f\n', R);
fprintf(fid, 'K_lqr = [%.8f %.8f %.8f %.8f]\n', K_lqr);
fprintf(fid, 'Nbar = %.8f\n', Nbar);
fprintf(fid, 'rise_time_x_s = %.6f\n', riseTimeX);
fprintf(fid, 'settling_time_x_s = %.6f\n', settlingTimeX);
fprintf(fid, 'settling_time_phi_s = %.6f\n', settlingTimePhi);
fprintf(fid, 'max_abs_phi_rad = %.6f\n', overshootPhiRad);
fprintf(fid, 'closed_loop_eigenvalues =\n');
for idx = 1:numel(closedLoopEigenvalues)
    fprintf(fid, '  %.8f%+.8fi\n', real(closedLoopEigenvalues(idx)), imag(closedLoopEigenvalues(idx)));
end
fclose(fid);

fprintf('Topic 5 GREEN: Tr_x=%.4f s, Ts_x=%.4f s, Ts_phi=%.4f s, max|phi|=%.4f rad\n', ...
    riseTimeX, settlingTimeX, settlingTimePhi, overshootPhiRad);

function create_topic5_simulink_model(projectRoot, A, B, K_lqr, Nbar)
    modelName = 'topic5_lqr_state_feedback';
    modelPath = fullfile(projectRoot, 'matlab', [modelName '.slx']);
    screenshotPath = fullfile(projectRoot, 'figures', 'topic5', 'simulation_structure.png');

    if bdIsLoaded(modelName)
        close_system(modelName, 0);
    end
    if exist(modelPath, 'file')
        load_system(modelPath);
        set_param(modelName, 'StopTime', '5', 'ScreenColor', 'white');
        assert_required_blocks(modelName, {'Reference', 'Nbar', 'Sum', 'Plant', 'Demux', 'K'});
        set_param([modelName '/Reference'], 'Time', '0', 'Before', '0', 'After', '0.2');
        set_param([modelName '/Nbar'], 'Gain', sprintf('%.12g', Nbar));
        set_param([modelName '/Plant'], ...
            'A', mat2str(A, 12), 'B', mat2str(B, 12), ...
            'C', 'eye(4)', 'D', 'zeros(4,1)');
        set_param([modelName '/K'], 'Gain', mat2str(K_lqr, 12), 'Multiplication', 'Matrix(K*u)');
        apply_topic5_visual_layout(modelName);
        style_simulink_model(modelName);
        assert_orthogonal_lines(modelName);
        set_param(modelName, 'ZoomFactor', 'FitSystem');
        save_system(modelName, modelPath);
        export_simulink_image(modelName, screenshotPath);
        close_system(modelName, 0);
        return;
    end
    new_system(modelName);
    open_system(modelName);
    set_param(modelName, 'Location', [80, 80, 1660, 760], 'StopTime', '5', 'ScreenColor', 'white');

    add_block('simulink/Sources/Step', [modelName '/Reference'], ...
        'Time', '0', 'Before', '0', 'After', '0.2', ...
        'Position', [45, 145, 135, 215]);
    add_block('simulink/Math Operations/Gain', [modelName '/Nbar'], ...
        'Gain', sprintf('%.12g', Nbar), ...
        'Position', [190, 135, 330, 225]);
    add_block('simulink/Math Operations/Sum', [modelName '/Sum'], ...
        'Inputs', '+-', 'IconShape', 'round', ...
        'Position', [390, 155, 420, 205]);
    add_block('simulink/Continuous/State-Space', [modelName '/Plant'], ...
        'A', mat2str(A, 12), 'B', mat2str(B, 12), ...
        'C', 'eye(4)', 'D', 'zeros(4,1)', ...
        'Position', [520, 135, 700, 225]);
    add_block('simulink/Signal Routing/Demux', [modelName '/Demux'], ...
        'Outputs', '4', 'Position', [790, 125, 800, 235]);
    add_block('simulink/Math Operations/Gain', [modelName '/K'], ...
        'Gain', mat2str(K_lqr, 12), 'Multiplication', 'Matrix(K*u)', ...
        'Position', [500, 330, 700, 420]);
    add_block('simulink/Sinks/Scope', [modelName '/Scope'], ...
        'NumInputPorts', '2', 'Position', [960, 125, 1040, 215]);
    add_block('simulink/Sinks/Out1', [modelName '/x'], ...
        'Position', [1060, 130, 1095, 160]);
    add_block('simulink/Sinks/Out1', [modelName '/Phi'], ...
        'Position', [1060, 205, 1095, 235]);

    add_line(modelName, 'Reference/1', 'Nbar/1', 'autorouting', 'on');
    add_line(modelName, 'Nbar/1', 'Sum/1', 'autorouting', 'on');
    add_line(modelName, 'Sum/1', 'Plant/1', 'autorouting', 'on');
    add_line(modelName, 'Plant/1', 'Demux/1', 'autorouting', 'on');
    add_line(modelName, 'Plant/1', 'K/1', 'autorouting', 'on');
    add_line(modelName, 'K/1', 'Sum/2', 'autorouting', 'on');
    add_line(modelName, 'Demux/1', 'Scope/1', 'autorouting', 'on');
    add_line(modelName, 'Demux/3', 'Scope/2', 'autorouting', 'on');
    add_line(modelName, 'Demux/1', 'x/1', 'autorouting', 'on');
    add_line(modelName, 'Demux/3', 'Phi/1', 'autorouting', 'on');

    apply_topic5_visual_layout(modelName);
    style_simulink_model(modelName);
    assert_orthogonal_lines(modelName);
    set_param(modelName, 'ZoomFactor', 'FitSystem');
    save_system(modelName, modelPath);
    export_simulink_image(modelName, screenshotPath);
    close_system(modelName, 0);
end

function export_simulink_image(modelName, screenshotPath)
    try
        Simulink.BlockDiagram.exportToImage(modelName, screenshotPath);
    catch
        print(['-s' modelName], '-dpng', '-r300', screenshotPath);
    end
    annotate_simulink_icons(modelName, screenshotPath);
    assert_white_image_background(screenshotPath);
end

function export_pole_figure(fig, outputPath, exportDpi)
    style_report_figure(fig);
    ax = findall(fig, 'Type', 'axes');
    lighten_sgrid(ax(1));
    poleLines = findall(ax(1), 'Type', 'line', 'Marker', 'x');
    if ~isempty(poleLines)
        set(poleLines, 'Color', [0.00, 0.22, 0.80], 'LineWidth', 2.8, 'MarkerSize', 12);
    end
    exportgraphics(fig, outputPath, 'Resolution', exportDpi, 'BackgroundColor', 'white', 'ContentType', 'image');
    assert_white_image_background(outputPath);
end

function annotate_simulink_icons(modelName, screenshotPath)
    rgb = imread(screenshotPath);
    if ndims(rgb) == 2
        rgb = repmat(rgb, 1, 1, 3);
    elseif size(rgb, 3) > 3
        rgb = rgb(:, :, 1:3);
    end
    [height, width, ~] = size(rgb);
    switch modelName
        case 'topic5_lqr_state_feedback'
            rgb = draw_step_icon(rgb, [0.028, 0.345, 0.060, 0.070], width, height);
            rgb = draw_text_icon(rgb, 'Nbar', [0.172, 0.318, 0.090, 0.065], width, height);
            rgb = draw_text_icon(rgb, '-K', [0.490, 0.846, 0.080, 0.055], width, height);
            rgb = draw_scope_icon(rgb, [0.922, 0.292, 0.050, 0.070], width, height);
    end
    imwrite(rgb, screenshotPath);
end

function rgb = draw_text_icon(rgb, label, box, width, height)
    pos = rel_box(box, width, height);
    rgb = insertText(rgb, [pos(1) + pos(3) / 2, pos(2) + pos(4) / 2], label, ...
        'AnchorPoint', 'Center', 'FontSize', max(22, round(pos(4) * 0.28)), ...
        'TextColor', 'black', 'BoxOpacity', 0);
end

function rgb = draw_step_icon(rgb, box, width, height)
    pos = rel_box(box, width, height);
    x = [pos(1) + 0.18 * pos(3), pos(1) + 0.45 * pos(3), pos(1) + 0.45 * pos(3), pos(1) + 0.82 * pos(3)];
    y = [pos(2) + 0.68 * pos(4), pos(2) + 0.68 * pos(4), pos(2) + 0.32 * pos(4), pos(2) + 0.32 * pos(4)];
    rgb = insertShape(rgb, 'Line', [x(1), y(1), x(2), y(2); x(2), y(2), x(3), y(3); x(3), y(3), x(4), y(4)], ...
        'Color', 'black', 'LineWidth', 3);
end

function rgb = draw_scope_icon(rgb, box, width, height)
    pos = rel_box(box, width, height);
    rgb = insertShape(rgb, 'Rectangle', [pos(1) + 0.18 * pos(3), pos(2) + 0.20 * pos(4), 0.64 * pos(3), 0.58 * pos(4)], ...
        'Color', 'black', 'LineWidth', 3);
    for k = 1:2
        x = pos(1) + (0.18 + 0.64 * k / 3) * pos(3);
        rgb = insertShape(rgb, 'Line', [x, pos(2) + 0.20 * pos(4), x, pos(2) + 0.78 * pos(4)], 'Color', [80, 80, 80], 'LineWidth', 1);
    end
    y = pos(2) + 0.49 * pos(4);
    rgb = insertShape(rgb, 'Line', [pos(1) + 0.18 * pos(3), y, pos(1) + 0.82 * pos(3), y], 'Color', [80, 80, 80], 'LineWidth', 1);
end

function pos = rel_box(box, width, height)
    pos = [box(1) * width, box(2) * height, box(3) * width, box(4) * height];
end

function lighten_sgrid(ax)
    gridLines = findall(ax, 'Type', 'line');
    for idx = 1:numel(gridLines)
        marker = get(gridLines(idx), 'Marker');
        if ~strcmp(marker, 'x')
            set(gridLines(idx), 'Color', [0.86, 0.86, 0.86], 'LineWidth', 0.65);
        end
    end
end

function apply_topic5_visual_layout(modelName)
    set_param([modelName '/Reference'], 'Position', [45, 145, 135, 215]);
    set_param([modelName '/Nbar'], 'Position', [190, 135, 330, 225]);
    set_param([modelName '/K'], 'Position', [500, 330, 700, 420]);
    set_param([modelName '/Scope'], 'Position', [960, 125, 1040, 215]);
    force_icon_display([modelName '/Reference']);
    force_icon_display([modelName '/Nbar']);
    force_icon_display([modelName '/K']);
    force_icon_display([modelName '/Scope']);
end

function force_icon_display(blockPath)
    params = get_param(blockPath, 'ObjectParameters');
    if isfield(params, 'IconDisplay')
        try
            set_param(blockPath, 'IconDisplay', 'Graphic');
        catch
        end
    end
    if isfield(params, 'ShowName')
        set_param(blockPath, 'ShowName', 'on');
    end
end

function configure_report_graphics()
    set(0, 'DefaultFigureVisible', 'off');
    set(0, 'DefaultFigureColor', 'w');
    set(0, 'DefaultAxesColor', 'w');
    set(0, 'DefaultAxesXColor', 'k');
    set(0, 'DefaultAxesYColor', 'k');
    set(0, 'DefaultAxesZColor', 'k');
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultAxesFontSize', 12);
    set(0, 'DefaultAxesLineWidth', 1.2);
    set(0, 'DefaultTextColor', 'k');
    set(0, 'DefaultTextFontName', 'Arial');
    set(0, 'DefaultTextFontSize', 12);
    set(0, 'DefaultLineLineWidth', 2.2);
end

function export_report_figure(fig, outputPath, exportDpi)
    style_report_figure(fig);
    exportgraphics(fig, outputPath, 'Resolution', exportDpi, 'BackgroundColor', 'white', 'ContentType', 'image');
    assert_white_image_background(outputPath);
end

function style_report_figure(fig)
    set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
    axesHandles = findall(fig, 'Type', 'axes');
    for idx = 1:numel(axesHandles)
        set(axesHandles(idx), 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
            'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
            'GridColor', [0.25, 0.25, 0.25], 'MinorGridColor', [0.55, 0.55, 0.55]);
    end
    textHandles = findall(fig, 'Type', 'text');
    if ~isempty(textHandles)
        set(textHandles, 'Color', 'k', 'FontName', 'Arial');
    end
    lineHandles = findall(fig, 'Type', 'line');
    if ~isempty(lineHandles)
        set(lineHandles, 'LineWidth', 2.2);
    end
    constantLines = findall(fig, 'Type', 'ConstantLine');
    if ~isempty(constantLines)
        set(constantLines, 'LineWidth', 1.9, 'Color', 'k');
    end
    legends = findall(fig, 'Type', 'legend');
    if ~isempty(legends)
        set(legends, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k', 'FontName', 'Arial', 'FontSize', 11);
    end
end

function style_simulink_model(modelName)
    set_param(modelName, 'ScreenColor', 'white');
    blockHandles = find_system(modelName, 'FindAll', 'on', 'Type', 'block');
    for idx = 1:numel(blockHandles)
        params = get_param(blockHandles(idx), 'ObjectParameters');
        if isfield(params, 'ForegroundColor')
            set_param(blockHandles(idx), 'ForegroundColor', 'black');
        end
        if isfield(params, 'BackgroundColor')
            set_param(blockHandles(idx), 'BackgroundColor', 'white');
        end
        if isfield(params, 'FontSize')
            set_param(blockHandles(idx), 'FontSize', '14');
        end
    end
end

function assert_required_blocks(modelName, blockNames)
    for idx = 1:numel(blockNames)
        blockPath = [modelName '/' blockNames{idx}];
        assert(getSimulinkBlockHandle(blockPath) ~= -1, ['Required block missing: ' blockPath]);
    end
end

function postprocess_simulink_image(imagePath)
    rgb = imread(imagePath);
    if ndims(rgb) == 2
        rgb = repmat(rgb, 1, 1, 3);
    elseif size(rgb, 3) > 3
        rgb = rgb(:, :, 1:3);
    end
    if ~isa(rgb, 'uint8')
        rgb = uint8(255 * double(rgb) / double(intmax(class(rgb))));
    end
    gray = 0.299 * double(rgb(:, :, 1)) + 0.587 * double(rgb(:, :, 2)) + 0.114 * double(rgb(:, :, 3));
    cornerSamples = [gray(1:20, 1:20); gray(1:20, end-19:end); gray(end-19:end, 1:20); gray(end-19:end, end-19:end)];
    output = uint8(255 * ones(size(rgb)));
    if median(cornerSamples(:)) < 80
        contentMask = gray > 100;
    else
        contentMask = gray < 190;
    end
    thickMask = conv2(double(contentMask), ones(5), 'same') > 0;
    for channelIdx = 1:3
        channel = output(:, :, channelIdx);
        channel(thickMask) = 0;
        output(:, :, channelIdx) = channel;
    end
    imwrite(output, imagePath);
end

function assert_white_image_background(imagePath)
    rgb = imread(imagePath);
    if ndims(rgb) == 2
        rgb = repmat(rgb, 1, 1, 3);
    elseif size(rgb, 3) > 3
        rgb = rgb(:, :, 1:3);
    end
    gray = 0.299 * double(rgb(:, :, 1)) + 0.587 * double(rgb(:, :, 2)) + 0.114 * double(rgb(:, :, 3));
    cornerSamples = [gray(1:20, 1:20); gray(1:20, end-19:end); gray(end-19:end, 1:20); gray(end-19:end, end-19:end)];
    assert(mean(cornerSamples(:)) > 245, 'Exported image background is not white.');
end

function assert_orthogonal_lines(modelName)
    lineHandles = find_system(modelName, 'FindAll', 'on', 'Type', 'line');
    for idx = 1:numel(lineHandles)
        points = get_param(lineHandles(idx), 'Points');
        for pointIdx = 1:(size(points, 1) - 1)
            dx = abs(points(pointIdx + 1, 1) - points(pointIdx, 1));
            dy = abs(points(pointIdx + 1, 2) - points(pointIdx, 2));
            assert(dx == 0 || dy == 0, 'Simulink line is not routed with right angles.');
        end
    end
end

function tr = rise_time(t, y, lowLevel, highLevel)
    lowIndex = find(y >= lowLevel, 1, 'first');
    highIndex = find(y >= highLevel, 1, 'first');
    if isempty(lowIndex) || isempty(highIndex)
        tr = inf;
    else
        tr = t(highIndex) - t(lowIndex);
    end
end

function ts = settling_time(t, y, target, tolerance)
    outside = find(abs(y - target) > tolerance);
    if isempty(outside)
        ts = 0;
    elseif outside(end) == numel(t)
        ts = inf;
    else
        ts = t(outside(end) + 1);
    end
end
