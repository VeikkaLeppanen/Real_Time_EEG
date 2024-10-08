#ifndef GACORRECTION_H
#define GACORRECTION_H

#include <mutex>
#include <cstddef>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

class GACorrection {
public:
    GACorrection(int    channel_count, 
                 int    template_avg_length,
                 int    TA_length)
                 :      channel_count_(channel_count),
                        template_avg_length_(template_avg_length),
                        TA_length_(TA_length),
                        number_of_samples_(template_avg_length * TA_length)
    {
        sample_data_ = Eigen::MatrixXd::Zero(channel_count_, number_of_samples_);
        correction_template_ = Eigen::MatrixXd::Zero(channel_count_, TA_length);
    }

    void update_template(int template_index, const Eigen::VectorXd &samples);
    void reset_index() { sample_data_index_ = 0; }

    Eigen::VectorXd getTemplateCol(int template_index) { return correction_template_.col(template_index); }
    Eigen::MatrixXd getTemplateCols(int template_index, int num_bundles);
    int getTemplateSize() { return correction_template_.cols(); }

    void printTemplate() { std::cout << correction_template_ << '\n'; }

private:
    Eigen::MatrixXd sample_data_;
    Eigen::MatrixXd correction_template_;
    size_t sample_data_index_ = 0;
    int channel_count_;
    int template_avg_length_;
    int TA_length_;
    int number_of_samples_;
};

#endif // GACORRECTION_H
