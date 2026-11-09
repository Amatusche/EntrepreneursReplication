"Empirical targets used by the paper calibration and model-fit tables."
function calibration_targets()
    _, DeBackerDataFull = makedebacker(false) # Exclude exits, which the model does not represent.
    DeBackerTargets = debacker_mobility_stats(DeBackerDataFull)
    Targets = CalibTargets(2.92, 0.26, 0.1, 0.2, 0.15, 0.46, 0.22, 0.22, 0.16, 7.1, 3.3, 6 / 10,
        DeBackerTargets.immobility_ratio, DeBackerTargets.rank_autocorrelation,
        DeBackerTargets.large_loss_probability,
        DeBackerTargets.large_gain_probability, DeBackerTargets.avg_deciles_lost) # KMV + negative liquid assets + entrepreneurial income mobility
    return Targets
end
