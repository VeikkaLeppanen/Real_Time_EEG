#include "GACorrection.h"


// TODO: multiple bundles
Eigen::MatrixXd GACorrection::getTemplateCols(int template_index, int num_bundles) {
    int cols_to_take = (template_index + num_bundles >= number_of_samples_) ? number_of_samples_ - 1 : template_index + num_bundles;

    return correction_template_.middleCols(template_index, cols_to_take); 
}

void GACorrection::update_template(int template_index, const Eigen::VectorXd &samples) {

    Eigen::VectorXd Old_avg = correction_template_.col(template_index);

    Eigen::VectorXd New_avg = Old_avg + (samples - sample_data_.col(sample_data_index_)) / template_avg_length_;

    correction_template_.col(template_index) = New_avg;

    sample_data_.col(sample_data_index_) = samples;
    sample_data_index_ = (sample_data_index_ + 1) % number_of_samples_;
    return; 
}