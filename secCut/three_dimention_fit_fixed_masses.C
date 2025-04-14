#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TFile.h"

void three_dimension_fit_fixed_masses() {
    // 定义变量
    RooRealVar m_Jpsi_1("m_Jpsi_1", "J/#psi_{1} invariant mass", 3.096, 3.0, 3.2);
    RooRealVar m_Jpsi_2("m_Jpsi_2", "J/#psi_{2} invariant mass", 3.096, 3.0, 3.2);
    RooRealVar m_Ups("m_Ups", "#Upsilon invariant mass", 9.460, 9.0, 10.0);

    // 定义数据集（假设你已经有一个包含数据的文件）
    RooDataSet data("data", "Dataset with m_Jpsi_1, m_Jpsi_2, m_Ups", RooArgSet(m_Jpsi_1, m_Jpsi_2, m_Ups));
    // 这里假设数据已经加载到 `data` 中，你可以从文件或手动填充数据

    // 定义信号和本底的PDF
    // J/psi_1 的高斯分布
    RooRealVar mean_Jpsi_1("mean_Jpsi_1", "Mean of J/#psi_{1}", 3.096); // 固定为 PDG 值
    mean_Jpsi_1.setConstant(true); // 固定不变质量的均值
    RooRealVar sigma_Jpsi_1("sigma_Jpsi_1", "Sigma of J/#psi_{1}", 0.01, 0.001, 0.1);
    RooGaussian gauss_Jpsi_1("gauss_Jpsi_1", "Gaussian for J/#psi_{1}", m_Jpsi_1, mean_Jpsi_1, sigma_Jpsi_1);

    // J/psi_2 的高斯分布
    RooRealVar mean_Jpsi_2("mean_Jpsi_2", "Mean of J/#psi_{2}", 3.096); // 固定为 PDG 值
    mean_Jpsi_2.setConstant(true); // 固定不变质量的均值
    RooRealVar sigma_Jpsi_2("sigma_Jpsi_2", "Sigma of J/#psi_{2}", 0.01, 0.001, 0.1);
    RooGaussian gauss_Jpsi_2("gauss_Jpsi_2", "Gaussian for J/#psi_{2}", m_Jpsi_2, mean_Jpsi_2, sigma_Jpsi_2);

    // Upsilon 的高斯分布
    RooRealVar mean_Ups("mean_Ups", "Mean of #Upsilon", 9.460); // 固定为 PDG 值
    mean_Ups.setConstant(true); // 固定不变质量的均值
    RooRealVar sigma_Ups("sigma_Ups", "Sigma of #Upsilon", 0.1, 0.01, 1.0);
    RooGaussian gauss_Ups("gauss_Ups", "Gaussian for #Upsilon", m_Ups, mean_Ups, sigma_Ups);

    // 定义联合信号PDF
    RooProdPdf signal_pdf("signal_pdf", "Signal PDF", RooArgSet(gauss_Jpsi_1, gauss_Jpsi_2, gauss_Ups));

    // 定义本底PDF（假设为均匀分布）
    RooUniform bkg_pdf("bkg_pdf", "Background PDF", RooArgSet(m_Jpsi_1, m_Jpsi_2, m_Ups));

    // 定义信号和本底的比例
    RooRealVar nsig("nsig", "Number of signal events", 1000, 0, 10000);
    RooRealVar nbkg("nbkg", "Number of background events", 1000, 0, 10000);

    // 定义总PDF
    RooAddPdf model("model", "Signal + Background", RooArgList(signal_pdf, bkg_pdf), RooArgList(nsig, nbkg));

    // 进行拟合
    RooFitResult* fitResult = model.fitTo(data, RooFit::Save());
    fitResult->Print(); // 打印拟合结果

    // 绘制结果
    TCanvas* c1 = new TCanvas("c1", "Fit Result: J/psi_1", 800, 600);
    RooPlot* frame_Jpsi_1 = m_Jpsi_1.frame();
    data.plotOn(frame_Jpsi_1);
    model.plotOn(frame_Jpsi_1);
    frame_Jpsi_1->Draw();

    TCanvas* c2 = new TCanvas("c2", "Fit Result: J/psi_2", 800, 600);
    RooPlot* frame_Jpsi_2 = m_Jpsi_2.frame();
    data.plotOn(frame_Jpsi_2);
    model.plotOn(frame_Jpsi_2);
    frame_Jpsi_2->Draw();

    TCanvas* c3 = new TCanvas("c3", "Fit Result: Upsilon", 800, 600);
    RooPlot* frame_Ups = m_Ups.frame();
    data.plotOn(frame_Ups);
    model.plotOn(frame_Ups);
    frame_Ups->Draw();

    // 保存结果
    c1->SaveAs("fit_result_Jpsi_1_fixed.png");
    c2->SaveAs("fit_result_Jpsi_2_fixed.png");
    c3->SaveAs("fit_result_Ups_fixed.png");
}