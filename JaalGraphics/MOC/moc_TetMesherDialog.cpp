/****************************************************************************
** Meta object code from reading C++ file 'TetMesherDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "TetMesherDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'TetMesherDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JTetMesherDialog_t {
    QByteArrayData data[10];
    char stringdata0[147];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JTetMesherDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JTetMesherDialog_t qt_meta_stringdata_JTetMesherDialog = {
    {
        QT_MOC_LITERAL(0, 0, 16), // "JTetMesherDialog"
        QT_MOC_LITERAL(1, 17, 10), // "genNewMesh"
        QT_MOC_LITERAL(2, 28, 0), // ""
        QT_MOC_LITERAL(3, 29, 11), // "fromHexMesh"
        QT_MOC_LITERAL(4, 41, 23), // "openTetGenOptionsDialog"
        QT_MOC_LITERAL(5, 65, 20), // "getBoundedDistortion"
        QT_MOC_LITERAL(6, 86, 14), // "getMesquiteOpt"
        QT_MOC_LITERAL(7, 101, 13), // "getStellarOpt"
        QT_MOC_LITERAL(8, 115, 19), // "openTetViewerDialog"
        QT_MOC_LITERAL(9, 135, 11) // "closeDialog"

    },
    "JTetMesherDialog\0genNewMesh\0\0fromHexMesh\0"
    "openTetGenOptionsDialog\0getBoundedDistortion\0"
    "getMesquiteOpt\0getStellarOpt\0"
    "openTetViewerDialog\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JTetMesherDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    8,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   54,    2, 0x08 /* Private */,
    3,    0,   55,    2, 0x08 /* Private */,
    4,    0,   56,    2, 0x08 /* Private */,
    5,    0,   57,    2, 0x08 /* Private */,
    6,    0,   58,    2, 0x08 /* Private */,
    7,    0,   59,    2, 0x08 /* Private */,
    8,    0,   60,    2, 0x08 /* Private */,
    9,    0,   61,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JTetMesherDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JTetMesherDialog *_t = static_cast<JTetMesherDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->genNewMesh();
            break;
        case 1:
            _t->fromHexMesh();
            break;
        case 2:
            _t->openTetGenOptionsDialog();
            break;
        case 3:
            _t->getBoundedDistortion();
            break;
        case 4:
            _t->getMesquiteOpt();
            break;
        case 5:
            _t->getStellarOpt();
            break;
        case 6:
            _t->openTetViewerDialog();
            break;
        case 7:
            _t->closeDialog();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JTetMesherDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JTetMesherDialog.data,
        qt_meta_data_JTetMesherDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JTetMesherDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JTetMesherDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JTetMesherDialog.stringdata0))
        return static_cast<void*>(const_cast< JTetMesherDialog*>(this));
    if (!strcmp(_clname, "Ui::TetMesherDialog"))
        return static_cast< Ui::TetMesherDialog*>(const_cast< JTetMesherDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JTetMesherDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
