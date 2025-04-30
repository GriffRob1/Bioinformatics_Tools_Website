import React from 'react';
import dna_code_image from '../images/dna-code-image.jpg'
import { useNavigate } from "react-router-dom";
import BlueButton from "./BlueButton";

const sectionStyle = {
    backgroundImage: 'url(' + dna_code_image + ')'
}

export default function TitleArea() {

    return (
        <section className={'container title-container'} style={sectionStyle}>
            <div className={'title-area'}>
                <h1 className={'title'}>Bioinformatics Tools</h1>
                <p>Essential tools for analysing biological data</p>
                <BlueButton content={'What is Bioinformatics?'}
                            URL={'/about'}
                            buttonClass={'what-is-bioinformatics-button'}/>
            </div>
        </section>
    )
}