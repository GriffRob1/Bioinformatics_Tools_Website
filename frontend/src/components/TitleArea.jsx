import React from 'react';
import dna_code_image from '../images/dna-code-image.jpg'

const sectionStyle = {
    backgroundImage: 'url(' + dna_code_image + ')',
}

export default function TitleArea() {
    return (
        <section className='container title-container' style={sectionStyle}>
            <div className='title-area'>
                <h1 className='title'>Bioinformatics Tools</h1>
            </div>
        </section>
    )
}